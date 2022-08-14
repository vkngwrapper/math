package math

import "golang.org/x/sys/cpu"

var supportsSSE bool

func init() {
	supportsSSE = cpu.X86.HasSSE2
}
