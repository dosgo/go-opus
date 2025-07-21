package opus

import (
	"errors"

	"github.com/dosgo/goOpus/celt"
	"github.com/dosgo/goOpus/comm"
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

	//gainDb := 0
	streams := 1
	coupledStreams := 0
	//var mapping []byte
	channelMap := false

	if len(d.Extradata) >= OPUS_HEAD_SIZE {
		//gainDb = int(binary.LittleEndian.Uint16(d.Extradata[16:18]))
		channelMap = d.Extradata[18] != 0
	}

	if len(d.Extradata) >= OPUS_HEAD_SIZE+2+channels {
		streams = int(d.Extradata[OPUS_HEAD_SIZE])
		coupledStreams = int(d.Extradata[OPUS_HEAD_SIZE+1])
		if streams+coupledStreams != channels {
			return errors.New("invalid channel mapping")
		}
		//mapping = d.Extradata[OPUS_HEAD_SIZE+2:]
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
	if opusPkt.Mode != comm.ModeCelt {
		d.Silk.Setup(opusPkt)
	}

	if opusPkt.Mode == comm.ModeCelt {
		d.Celt.Setup(opusPkt.Stereo)
	}

	if opusPkt.Mode == comm.ModeHybrid {
		// TODO: 实现混合模式
	}

	// 解码所有帧
	for _, frame := range opusPkt.Frames {
		rd := comm.NewRangeDecoder(frame)

		if opusPkt.Mode != comm.ModeCelt {
			if _, err := d.Silk.Decode(rd); err != nil {
				return err
			}
		} else {
			d.Silk.Flush()
		}

		size := len(frame)
		consumed := rd.Tell()
		redundancy := false

		if opusPkt.Mode == comm.ModeHybrid && consumed+37 <= size*8 {
			redundancy = rd.DecodeLogP(12)
		} else if opusPkt.Mode == comm.ModeSilk && consumed+17 <= size*8 {
			redundancy = true
		}

		if redundancy {
			redundancyPos := rd.DecodeLogP(1)
			redundancySize := 0

			if opusPkt.Mode == comm.ModeHybrid {
				redundancySize = int(rd.DecodeUniform(256)) + 2
			} else {
				redundancySize = size - (consumed+7)/8
			}

			if redundancySize >= size {
				return errors.New("invalid redundancy size")
			}

			//_size := size - redundancySize

			if redundancyPos {
				// TODO: 实现冗余解码
				//d.Celt.Flush()
			}
		}

		if opusPkt.Mode != comm.ModeSilk {
			outBuf := make([]float32, 1024)
			var bandStart = 0
			var bandEnd = 0
			if opusPkt.Mode == comm.ModeHybrid {

				bandStart = ternary(opusPkt.Mode == comm.ModeHybrid, 17, 0)
				bandEnd = opusPkt.Bandwidth.CeltBand()

			}
			d.Celt.Decode(rd, outBuf, opusPkt.FrameDuration, bandStart, bandEnd)
		}
	}

	return nil
}
func ternary(condition bool, trueVal, falseVal int) int {
	if condition {
		return trueVal
	}
	return falseVal
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
		//d.Celt.Flush()
	}
	return nil
}

// ParsePacket 解析Opus包
func ParsePacket(data []byte) (*comm.Packet, error) {
	// 简化的包解析实现
	pkt := &comm.Packet{
		FrameDuration: comm.FrameDurationStandard,
		Stereo:        len(data) > 100, // 简化判断
		Bandwidth:     comm.BandwidthWide,
		Mode:          comm.ModeHybrid,
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
