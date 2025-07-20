package opus

import (
	"encoding/binary"
	"errors"
	"fmt"
	"os"

	"github.com/dosgo/goOpus/celt"
	"github.com/dosgo/goOpus/silk"
)

const OPUS_HEAD_SIZE = 19

// Decoder 表示 Opus 音频解码器
type Decoder struct {
	Extradata []byte
	Silk      *silk.Silk
	Celt      *celt.Celt
}

// NewDecoder 创建新的 Opus 解码器
func NewDecoder() *Decoder {
	return &Decoder{}
}

// SetExtradata 设置额外的配置数据
func (d *Decoder) SetExtradata(extra []byte) {
	d.Extradata = make([]byte, len(extra))
	copy(d.Extradata, extra)
}

// Configure 配置解码器
func (d *Decoder) Configure() error {
	if len(d.Extradata) == 0 {
		return errors.New("missing extradata")
	}

	channels := 2
	if len(d.Extradata) > 9 {
		channels = int(d.Extradata[9])
	}

	gainDb := 0
	streams := 1
	coupledStreams := 0
	var mapping []byte
	channelMap := false

	if len(d.Extradata) >= OPUS_HEAD_SIZE {
		gainDb = int(binary.LittleEndian.Uint16(d.Extradata[16:18]))
		channelMap = d.Extradata[18] != 0
	}

	if len(d.Extradata) >= OPUS_HEAD_SIZE+2+channels {
		streams = int(d.Extradata[OPUS_HEAD_SIZE])
		coupledStreams = int(d.Extradata[OPUS_HEAD_SIZE+1])
		if streams+coupledStreams != channels {
			return errors.New("invalid channel mapping")
		}
		mapping = d.Extradata[OPUS_HEAD_SIZE+2:]
	} else {
		if channels > 2 || channelMap {
			return errors.New("configuration invalid")
		}
		if channels > 1 {
			coupledStreams = 1
		}
	}

	stereo := channels > 1
	d.Silk = silk.NewSilk(stereo)
	d.Celt = celt.NewCelt(stereo)

	return nil
}

// SendPacket 发送音频包进行解码
func (d *Decoder) SendPacket(packetData []byte) error {
	if d.Silk == nil || d.Celt == nil {
		return errors.New("decoder not configured")
	}

	opusPkt, err := ParsePacket(packetData)
	if err != nil {
		return err
	}

	// 配置 SILK 和 CELT 解码器
	if opusPkt.Mode != ModeCelt {
		d.Silk.Setup(opusPkt)
	}

	if opusPkt.Mode == ModeCelt {
		d.Celt.Setup(opusPkt)
	}

	if opusPkt.Mode == ModeHybrid {
		// TODO: 实现混合模式
	}

	// 解码所有帧
	for _, frame := range opusPkt.Frames {
		rd := NewRangeDecoder(frame)

		if opusPkt.Mode != ModeCelt {
			if err := d.Silk.Decode(rd); err != nil {
				return err
			}
		} else {
			d.Silk.Flush()
		}

		size := len(frame)
		consumed := rd.Tell()
		redundancy := false

		if opusPkt.Mode == ModeHybrid && consumed+37 <= size*8 {
			redundancy = rd.DecodeLogP(12) != 0
		} else if opusPkt.Mode == ModeSilk && consumed+17 <= size*8 {
			redundancy = true
		}

		if redundancy {
			redundancyPos := rd.DecodeLogP(1) != 0
			redundancySize := 0

			if opusPkt.Mode == ModeHybrid {
				redundancySize = int(rd.DecodeUniform(256)) + 2
			} else {
				redundancySize = size - (consumed+7)/8
			}

			if redundancySize >= size {
				return errors.New("invalid redundancy size")
			}

			_size := size - redundancySize

			if redundancyPos {
				// TODO: 实现冗余解码
				d.Celt.Flush()
			}
		}

		if opusPkt.Mode != ModeSilk {
			outBuf := make([]float32, 1024)
			bandRange := NewRange(0, 0)

			if opusPkt.Mode == ModeHybrid {
				bandRange = NewRange(17, opusPkt.Bandwidth.CeltBand())
			}

			if err := d.Celt.Decode(rd, outBuf, opusPkt.FrameDuration, bandRange); err != nil {
				return err
			}
		}
	}

	return nil
}

// ReceiveFrame 接收解码后的音频帧
func (d *Decoder) ReceiveFrame() ([]float32, []float32, error) {
	if d.Silk == nil {
		return nil, nil, errors.New("decoder not configured")
	}

	// 获取左声道数据
	left := make([]float32, len(d.Silk.LeftOutbuf))
	copy(left, d.Silk.LeftOutbuf)

	// 获取右声道数据
	right := make([]float32, len(d.Silk.RightOutbuf))
	copy(right, d.Silk.RightOutbuf)

	return left, right, nil
}

// Flush 刷新解码器状态
func (d *Decoder) Flush() error {
	if d.Silk != nil {
		d.Silk.Flush()
	}
	if d.Celt != nil {
		d.Celt.Flush()
	}
	return nil
}

// Test 测试函数
func Test(testVectors []string) {
	for i, file := range testVectors {
		testSendPacket(i+1, file)
	}
}

func testSendPacket(index int, filename string) {
	fmt.Printf("Testing %s\n", filename)

	// 读取整个文件内容
	data, err := os.ReadFile(filename)
	if err != nil {
		fmt.Printf("Error reading file: %v\n", err)
		return
	}

	// 创建解码器
	dec := NewDecoder()

	// 假设文件前19字节是extradata
	if len(data) >= OPUS_HEAD_SIZE {
		dec.SetExtradata(data[:OPUS_HEAD_SIZE])
	}

	// 配置解码器
	if err := dec.Configure(); err != nil {
		fmt.Printf("Configuration error: %v\n", err)
		return
	}

	// 发送包进行解码
	if err := dec.SendPacket(data); err != nil {
		fmt.Printf("Send packet error: %v\n", err)
	}

	// 获取解码后的帧
	left, right, err := dec.ReceiveFrame()
	if err != nil {
		fmt.Printf("Receive frame error: %v\n", err)
		return
	}

	fmt.Printf("Decoded %d left samples, %d right samples\n", len(left), len(right))
}

// ParsePacket 解析Opus包
func ParsePacket(data []byte) (*Packet, error) {
	// 简化的包解析实现
	pkt := &Packet{
		FrameDuration: FrameDurationStandard,
		Stereo:        len(data) > 100, // 简化判断
		Bandwidth:     BandwidthWide,
		Mode:          ModeHybrid,
	}

	// 假设包包含多个帧
	frameSize := len(data) / 3
	if frameSize > 0 {
		pkt.Frames = [][]byte{
			data[:frameSize],
			data[frameSize : 2*frameSize],
			data[2*frameSize:],
		}
	} else {
		pkt.Frames = [][]byte{data}
	}

	return pkt, nil
}

// Packet 表示Opus音频包
type Packet struct {
	FrameDuration FrameDuration
	Stereo        bool
	Bandwidth     Bandwidth
	Mode          Mode
	Frames        [][]byte
}

// FrameDuration 表示帧持续时间类型
type FrameDuration int

const (
	FrameDurationMedium FrameDuration = iota
	FrameDurationStandard
	FrameDurationLong
	FrameDurationVeryLong
)

// Bandwidth 表示带宽类型
type Bandwidth int

const (
	BandwidthNarrow Bandwidth = iota
	BandwidthMedium
	BandwidthWide
	BandwidthSuperwide
	BandwidthFull
)

// CeltBand 获取CELT带宽
func (b Bandwidth) CeltBand() int {
	switch b {
	case BandwidthNarrow:
		return 64
	case BandwidthMedium:
		return 96
	default:
		return 128
	}
}

// Mode 表示编码模式
type Mode int

const (
	ModeSilk Mode = iota
	ModeCelt
	ModeHybrid
)
