package comm

import (
	"errors"
)

// Code 表示数据包编码类型
type Code int

const (
	CodeSingle      Code = iota // 单帧数据包
	CodeDoubleEqual             // 等长双帧数据包
	CodeDoubleVary              // 变长双帧数据包
	CodeMultiple                // 多帧数据包
)

// Mode 表示编码模式
type Mode int

const (
	ModeSilk   Mode = iota // SILK 模式
	ModeCelt               // CELT 模式
	ModeHybrid             // 混合模式
)

// Bandwidth 表示带宽类型
type Bandwidth int

const (
	BandwidthNarrow    Bandwidth = 8000  // 窄带
	BandwidthMedium    Bandwidth = 12000 // 中带
	BandwidthWide      Bandwidth = 16000 // 宽带
	BandwidthSuperWide Bandwidth = 24000 // 超宽带
	BandwidthFull      Bandwidth = 48000 // 全带
)

// CeltBand 获取CELT带宽
func (b Bandwidth) CeltBand() int {
	switch b {
	case BandwidthNarrow:
		return 13
	case BandwidthMedium:
		return 17
	case BandwidthWide:
		return 17
	case BandwidthSuperWide:
		return 19
	case BandwidthFull:
		return 21
	default:
		return 17
	}
}

// FrameDuration 表示帧持续时间
type FrameDuration int

const (
	FrameDurationVeryShort FrameDuration = 120  // 2.5ms
	FrameDurationShort     FrameDuration = 240  // 5ms
	FrameDurationMedium    FrameDuration = 480  // 10ms
	FrameDurationStandard  FrameDuration = 960  // 20ms
	FrameDurationLong      FrameDuration = 1920 // 40ms
	FrameDurationVeryLong  FrameDuration = 2880 // 60ms
)

// Packet 表示Opus数据包
type Packet struct {
	Code          Code
	VBR           bool
	Config        int
	Stereo        bool
	Padding       int
	Mode          Mode
	Bandwidth     Bandwidth
	FrameDuration FrameDuration
	Frames        [][]byte
}

const (
	maxFrameSize = 1275
	maxFrames    = 48
	maxPacketDur = 5760
)

// ParsePacket 解析Opus数据包
func ParsePacket(buf []byte) (*Packet, error) {
	if len(buf) < 1 {
		return nil, errors.New("empty packet")
	}

	p := &Packet{
		Frames: make([][]byte, 0),
	}

	// 解析头部字节
	code := buf[0] & 0x03
	p.Config = int(buf[0]>>3) & 0x1F
	p.Stereo = (buf[0]>>2)&0x01 == 1

	// 移动到数据部分
	data := buf[1:]

	switch code {
	case 0: // 单帧数据包
		if err := p.parseSinglePacket(data); err != nil {
			return nil, err
		}
	case 1: // 等长双帧数据包
		if err := p.parseDoubleEqualPacket(data); err != nil {
			return nil, err
		}
	case 2: // 变长双帧数据包
		if err := p.parseDoubleVaryPacket(data); err != nil {
			return nil, err
		}
	case 3: // 多帧数据包
		if err := p.parseMultiplePacket(data); err != nil {
			return nil, err
		}
	default:
		return nil, errors.New("invalid packet code")
	}

	// 根据配置设置模式、带宽和帧持续时间
	if err := p.setConfig(); err != nil {
		return nil, err
	}

	return p, nil
}

// parseSinglePacket 解析单帧数据包
func (p *Packet) parseSinglePacket(data []byte) error {
	p.Code = CodeSingle
	p.VBR = false
	p.Frames = append(p.Frames, data)
	return nil
}

// parseDoubleEqualPacket 解析等长双帧数据包
func (p *Packet) parseDoubleEqualPacket(data []byte) error {
	p.Code = CodeDoubleEqual
	p.VBR = false

	if len(data)%2 != 0 {
		return errors.New("invalid double equal packet length")
	}

	half := len(data) / 2
	p.Frames = append(p.Frames, data[:half], data[half:])
	return nil
}

// parseDoubleVaryPacket 解析变长双帧数据包
func (p *Packet) parseDoubleVaryPacket(data []byte) error {
	p.Code = CodeDoubleVary
	p.VBR = true

	// 解析第一帧长度
	offset, len1, err := xiphLacingU16(data)
	if err != nil {
		return err
	}

	if len1+offset > len(data) {
		return errors.New("invalid double vary packet length")
	}

	// 分割两帧
	frame1 := data[offset : offset+len1]
	frame2 := data[offset+len1:]

	p.Frames = append(p.Frames, frame1, frame2)
	return nil
}

// parseMultiplePacket 解析多帧数据包
func (p *Packet) parseMultiplePacket(data []byte) error {
	p.Code = CodeMultiple
	p.VBR = (data[0]>>7)&0x01 == 1

	count := int(data[0] & 0x3F)
	padding := (data[0]>>6)&0x01 == 1

	if count == 0 || count > maxFrames {
		return errors.New("invalid frame count")
	}

	// 处理填充
	offset := 1
	if padding {
		off, pad, err := xiphLacingU32(data[1:])
		if err != nil {
			return err
		}
		p.Padding = pad
		offset += off
		data = data[:len(data)-pad]
	}

	// 移动到帧数据
	frameData := data[offset:]

	if p.VBR {
		// 变长帧
		lens := make([]int, count-1)
		for i := 0; i < count-1; i++ {
			off, len, err := xiphLacingU16(frameData)
			if err != nil {
				return err
			}
			lens[i] = len
			frameData = frameData[off:]
		}

		// 分割帧
		for _, _len := range lens {
			if _len > len(frameData) {
				return errors.New("invalid frame length")
			}
			p.Frames = append(p.Frames, frameData[:_len])
			frameData = frameData[_len:]
		}
		p.Frames = append(p.Frames, frameData) // 最后一帧
	} else {
		// 等长帧
		frameLen := len(frameData) / count
		if frameLen*count != len(frameData) || frameLen > maxFrameSize {
			return errors.New("invalid frame length")
		}

		for i := 0; i < count; i++ {
			start := i * frameLen
			end := start + frameLen
			p.Frames = append(p.Frames, frameData[start:end])
		}
	}

	return nil
}

// setConfig 根据配置设置模式、带宽和帧持续时间
func (p *Packet) setConfig() error {
	switch {
	case p.Config >= 0 && p.Config <= 11:
		p.Mode = ModeSilk
		switch {
		case p.Config <= 3:
			p.Bandwidth = BandwidthNarrow
		case p.Config <= 7:
			p.Bandwidth = BandwidthMedium
		case p.Config <= 11:
			p.Bandwidth = BandwidthWide
		}
		switch p.Config & 0x03 {
		case 0:
			p.FrameDuration = FrameDurationMedium
		case 1:
			p.FrameDuration = FrameDurationStandard
		case 2:
			p.FrameDuration = FrameDurationLong
		case 3:
			p.FrameDuration = FrameDurationVeryLong
		}
	case p.Config >= 12 && p.Config <= 15:
		p.Mode = ModeHybrid
		switch {
		case p.Config <= 13:
			p.Bandwidth = BandwidthSuperWide
		case p.Config <= 15:
			p.Bandwidth = BandwidthFull
		}
		switch p.Config & 0x01 {
		case 0:
			p.FrameDuration = FrameDurationMedium
		case 1:
			p.FrameDuration = FrameDurationStandard
		}
	case p.Config >= 16 && p.Config <= 31:
		p.Mode = ModeCelt
		switch {
		case p.Config <= 19:
			p.Bandwidth = BandwidthNarrow
		case p.Config <= 23:
			p.Bandwidth = BandwidthWide
		case p.Config <= 27:
			p.Bandwidth = BandwidthSuperWide
		case p.Config <= 31:
			p.Bandwidth = BandwidthFull
		}
		switch p.Config & 0x03 {
		case 0:
			p.FrameDuration = FrameDurationVeryShort
		case 1:
			p.FrameDuration = FrameDurationShort
		case 2:
			p.FrameDuration = FrameDurationMedium
		case 3:
			p.FrameDuration = FrameDurationStandard
		}
	default:
		return errors.New("invalid config value")
	}
	return nil
}

// xiphLacingU16 解析Xiph风格的16位长度编码
func xiphLacingU16(buf []byte) (offset int, length int, err error) {
	if len(buf) < 1 {
		return 0, 0, errors.New("buffer too short")
	}

	v := int(buf[0])
	if v < 252 {
		return 1, v, nil
	}

	if len(buf) < 2 {
		return 0, 0, errors.New("buffer too short for extended length")
	}

	length = v + 4*int(buf[1])
	return 2, length, nil
}

// xiphLacingU32 解析Xiph风格的32位长度编码
func xiphLacingU32(buf []byte) (offset int, length int, err error) {
	v := 0
	for i, b := range buf {
		v += int(b)
		offset = i + 1
		if b < 255 {
			break
		}
		v -= 1
	}
	return offset, v, nil
}
