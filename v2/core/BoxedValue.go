package opus

// BoxedValueByte 封装 byte 值以实现引用传递
type BoxedValueByte struct {
	Val byte
}

// NewBoxedValueByte 创建 BoxedValueByte 实例
func NewBoxedValueByte(v byte) *BoxedValueByte {
	return &BoxedValueByte{Val: v}
}

// BoxedValueShort 封装 int16 值以实现引用传递
type BoxedValueShort struct {
	Val int16
}

// NewBoxedValueShort 创建 BoxedValueShort 实例
func NewBoxedValueShort(v int16) *BoxedValueShort {
	return &BoxedValueShort{Val: v}
}

// BoxedValueInt 封装 int 值以实现引用传递
type BoxedValueInt struct {
	Val int
}

// NewBoxedValueInt 创建 BoxedValueInt 实例
func NewBoxedValueInt(v int) *BoxedValueInt {
	return &BoxedValueInt{Val: v}
}
