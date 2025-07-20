package opus

import (
	"encoding/binary"
	"math"
	"unsafe"
)

// 切片操作工具包

// Initialize2DInt 创建二维int切片
func Initialize2DInt(x, y int) [][]int {
	matrix := make([][]int, x)
	for i := range matrix {
		matrix[i] = make([]int, y)
	}
	return matrix
}

// Initialize2DFloat 创建二维float64切片
func Initialize2DFloat(x, y int) [][]float64 {
	matrix := make([][]float64, x)
	for i := range matrix {
		matrix[i] = make([]float64, y)
	}
	return matrix
}

// Initialize2DInt16 创建二维int16切片
func Initialize2DInt16(x, y int) [][]int16 {
	matrix := make([][]int16, x)
	for i := range matrix {
		matrix[i] = make([]int16, y)
	}
	return matrix
}

// Initialize2DByte 创建二维byte切片
func Initialize2DByte(x, y int) [][]byte {
	matrix := make([][]byte, x)
	for i := range matrix {
		matrix[i] = make([]byte, y)
	}
	return matrix
}

// Initialize3DByte 创建三维byte切片
func Initialize3DByte(x, y, z int) [][][]byte {
	matrix := make([][][]byte, x)
	for i := range matrix {
		matrix[i] = make([][]byte, y)
		for j := range matrix[i] {
			matrix[i][j] = make([]byte, z)
		}
	}
	return matrix
}

// SetByteSlice 填充整个byte切片
func SetByteSlice(slice []byte, value byte) {
	for i := range slice {
		slice[i] = value
	}
}

// SetInt16Slice 填充整个int16切片
func SetInt16Slice(slice []int16, value int16) {
	for i := range slice {
		slice[i] = value
	}
}

// SetIntSlice 填充整个int切片
func SetIntSlice(slice []int, value int) {
	for i := range slice {
		slice[i] = value
	}
}

// SetFloat32Slice 填充整个float32切片
func SetFloat32Slice(slice []float32, value float32) {
	for i := range slice {
		slice[i] = value
	}
}

// SetByteSliceLen 填充byte切片的前length个元素
func SetByteSliceLen(slice []byte, value byte, length int) {
	for i := 0; i < length && i < len(slice); i++ {
		slice[i] = value
	}
}

// SetInt16SliceLen 填充int16切片的前length个元素
func SetInt16SliceLen(slice []int16, value int16, length int) {
	for i := 0; i < length && i < len(slice); i++ {
		slice[i] = value
	}
}

// SetIntSliceLen 填充int切片的前length个元素
func SetIntSliceLen(slice []int, value int, length int) {
	for i := 0; i < length && i < len(slice); i++ {
		slice[i] = value
	}
}

// SetFloat32SliceLen 填充float32切片的前length个元素
func SetFloat32SliceLen(slice []float32, value float32, length int) {
	for i := 0; i < length && i < len(slice); i++ {
		slice[i] = value
	}
}

// SetByteSliceOffset 在byte切片的指定偏移和长度内填充值
func SetByteSliceOffset(slice []byte, value byte, offset, length int) {
	end := offset + length
	if end > len(slice) {
		end = len(slice)
	}

	for i := offset; i < end; i++ {
		slice[i] = value
	}
}

// SetInt16SliceOffset 在int16切片的指定偏移和长度内填充值
func SetInt16SliceOffset(slice []int16, value int16, offset, length int) {
	end := offset + length
	if end > len(slice) {
		end = len(slice)
	}

	for i := offset; i < end; i++ {
		slice[i] = value
	}
}

// SetIntSliceOffset 在int切片的指定偏移和长度内填充值
func SetIntSliceOffset(slice []int, value int, offset, length int) {
	end := offset + length
	if end > len(slice) {
		end = len(slice)
	}

	for i := offset; i < end; i++ {
		slice[i] = value
	}
}

// MoveByteSlice 在同一个byte切片内移动数据块
func MoveByteSlice(slice []byte, srcIdx, dstIdx, length int) {
	// 使用copy处理源和目标重叠的情况
	copy(slice[dstIdx:], slice[srcIdx:srcIdx+length])
}

// MoveInt16Slice 在同一个int16切片内移动数据块
func MoveInt16Slice(slice []int16, srcIdx, dstIdx, length int) {
	// 使用copy处理源和目标重叠的情况
	copy(slice[dstIdx:], slice[srcIdx:srcIdx+length])
}

// MoveIntSlice 在同一个int切片内移动数据块
func MoveIntSlice(slice []int, srcIdx, dstIdx, length int) {
	// 使用copy处理源和目标重叠的情况
	copy(slice[dstIdx:], slice[srcIdx:srcIdx+length])
}

// ======== 高效版本 (使用内存块操作) ========

// FastSetUint32Slice 高效填充uint32切片
func FastSetUint32Slice(slice []uint32, value uint32) {
	// 使用内存复制加速
	bytes := make([]byte, 4)
	binary.LittleEndian.PutUint32(bytes, value)

	block := []byte{bytes[0], bytes[1], bytes[2], bytes[3]}
	for i := 0; i < len(slice); i += 4 {
		start := i * 4
		end := start + 4
		if end > len(slice)*4 {
			end = len(slice) * 4
		}
		copy(slice[i:], bytes)
	}
}

// FastSetFloat32Slice 高效填充float32切片
func FastSetFloat32Slice(slice []float32, value float32) {
	// 转换为uint32再使用内存复制
	bits := math.Float32bits(value)
	FastSetUint32Slice(
		*(*[]uint32)(unsafe.Pointer(&slice)),
		bits,
	)
}

// 注：使用unsafe的方法需要导入unsafe包，并且只在性能敏感场景使用
