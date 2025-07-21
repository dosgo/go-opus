package main

import (
	"fmt"

	"github.com/dosgo/goOpus/opus"
)

func main() {
	ddd := opus.NewDecoder()
	for {
		ddd.SendPacket([]byte{})
		_, _, err := ddd.ReceiveFrame()
		if err != nil {
			fmt.Printf("ReceiveFrame err:%v\r\n", err)
		}

	}

}
