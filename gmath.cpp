#include <cstdint>

uint32_t mcn_mulm(uint32_t a, uint32_t b, uint32_t mod)
{
	uint32_t r;
	__asm__ __volatile__ (
		"movl %%ecx, %%eax\n\t"
		"mull %%edx\n\t"
		"divl %%r8d\n\t"
		"movl %%edx, %%eax\n\t"
		: "=a" (r)
	);
	return r;
}

uint64_t mcn_mulm(uint64_t a, uint64_t b, uint64_t mod)
{
	uint64_t r;
	__asm__ __volatile__ (
		"movq %%rcx, %%rax\n\t"
		"mulq %%rdx\n\t"
		"divq %%r8\n\t"
		"movq %%rdx, %%rax\n\t"
		: "=a" (r)
	);
	return r;
}

uint32_t mcn_addm(uint32_t a, uint32_t b, uint32_t mod)
{
	uint32_t r;
	__asm__ __volatile__ (
		"movl %%ecx, %%eax\n\t"
		"addl %%edx, %%eax\n\t"
		"movl $0,    %%edx\n\t"
		"jnc  L1%=\n\t"
		"incl %%edx\n\t"
		"L1%=:\n\t"
		"divl %%r8d\n\t"
		"movl %%edx, %%eax\n\t"
		: "=a" (r)
	);
	return r;
}

uint32_t mcn_subm(uint32_t a, uint32_t b, uint32_t mod)
{
	uint32_t r;
	__asm__ __volatile__ (
		"movl %%ecx, %%eax\n\t"
		"subl %%edx, %%eax\n\t"
		"movl $0,    %%edx\n\t"
		"jnb  L1%=\n\t"
		"addl %%r8d, %%eax\n\t"
		"L1%=:\n\t"
		"divl %%r8d\n\t"
		"movl %%edx, %%eax\n\t"
		: "=a" (r)
	);
	return r;
}

uint64_t mcn_addm(uint64_t a, uint64_t b, uint64_t mod)
{
	uint64_t r;
	__asm__ __volatile__ (
		"movq %%rcx, %%rax\n\t"
		"addq %%rdx, %%rax\n\t"
		"movq $0,    %%rdx\n\t"
		"jnc  L1%=\n\t"
		"incq %%rdx\n\t"
		"L1%=:\n\t"
		"divq %%r8\n\t"
		"movq %%rdx, %%rax\n\t"
		: "=a" (r)
	);
	return r;
}

uint64_t mcn_subm(uint64_t a, uint64_t b, uint64_t mod)
{
	uint64_t r;
	__asm__ __volatile__ (
		"movq %%rcx, %%rax\n\t"
		"subq %%rdx, %%rax\n\t"
		"movq $0,    %%rdx\n\t"
		"jnb  L1%=\n\t"
		"addq %%r8,  %%rax\n\t"
		"L1%=:\n\t"
		"divq %%r8\n\t"
		"movq %%rdx, %%rax\n\t"
		: "=a" (r)
	);
	return r;
}

uint32_t mpn_mul_1(uint32_t *f, const uint32_t *g, int n, uint32_t h)
{
	uint32_t p;
	__asm__ __volatile__ (
		"xorl  %%ebx, %%ebx\n\t"
		"movl  %%edx, %%r8d\n\t"
		"L1%=:\n\t"
		"lodsl\n\t"
		"mull  %%r8d\n\t"
		"addl  %%ebx, %%eax\n\t"
		"adcl  $0, %%edx\n\t"
		"stosl\n\t"
		"movl  %%edx, %%ebx\n\t"
		"loop  L1%=\n\t"
		: "=b" (p)
		: "c" (n), "d" (h), "D" (f), "S" (g)
	);
	return p;
}

uint64_t mpn_mul_1(uint64_t *f, const uint64_t *g, long n, uint64_t h)
{
	uint64_t p;
	__asm__ __volatile__ (
		"xorq  %%rbx, %%rbx\n\t"
		"movq  %%rdx, %%r8\n\t"
		"cld\n\t"
	"L1%=:\n\t"
		"lodsq\n\t"
		"mulq  %%r8\n\t"
		"addq  %%rbx, %%rax\n\t"
		"adcq  $0, %%rdx\n\t"
		"stosq\n\t"
		"movq  %%rdx, %%rbx\n\t"
		"loop  L1%=\n\t"
		: "=b" (p)
		: "c" (n), "d" (h), "D" (f), "S" (g)
	);
	return p;
}

uint32_t mpn_add_1(uint32_t *f, const uint32_t *g, int n, uint32_t m)
{
	uint32_t p;
	__asm__ __volatile__ (
		"xorq  %%rbx, %%rbx\n\t"
		"L1%=:\n\t"
		"adcl (%%rsi, %%rbx, 4), %%eax\n\t"
		"movl  %%eax,(%%rdi, %%rbx, 4)\n\t"
		"movl  $0, %%eax\n\t"
		"incq  %%rbx\n\t"
		"loop  L1%=\n\t"
		"movl  $1, %%eax\n\t"
		"jc    L2%=\n\t"
		"xorl  %%eax, %%eax\n\t"
		"L2%=:\n\t"
		: "=a" (p) : "0" (m), "c" (n), "D" (f), "S" (g) : "%rbx"
	);
	return p;
}

uint32_t mpn_add_n (uint32_t *f, const uint32_t *g, const uint32_t *h, int n)
{
	uint32_t p;
	__asm__ __volatile__ (
		"xorq  %%rbx, %%rbx\n\t"
		"L1%=:\n\t"
		"movl (%%rsi, %%rbx, 4), %%eax\n\t"
		"adcl (%%rdx, %%rbx, 4), %%eax\n\t"
		"movl  %%eax,(%%rdi, %%rbx, 4)\n\t"
		"incq  %%rbx\n\t"
		"loop  L1%=\n\t"
		"movl  $1, %%eax\n\t"
		"jc    L2%=\n\t"
		"xorl  %%eax, %%eax\n\t"
		"L2%=:\n\t"
		: "=a" (p) : "c" (n), "d" (h), "S" (g), "D" (f) : "%rbx"
	);
	return p;
}

uint64_t mpn_add_1(uint64_t *f, const uint64_t *g, long n, uint64_t m)
{
	uint64_t p;
	__asm__ __volatile__ (
		"xorq  %%rbx, %%rbx\n\t"
		"L1%=:\n\t"
		"adcq (%%rsi, %%rbx, 8), %%rax\n\t"
		"movq  %%rax,(%%rdi, %%rbx, 8)\n\t"
		"movq  $0, %%rax\n\t"
		"incq  %%rbx\n\t"
		"loop  L1%=\n\t"
		"movq  $1, %%rax\n\t"
		"jc    L2%=\n\t"
		"xorq  %%rax, %%rax\n\t"
		"L2%=:\n\t"
		: "=a" (p) : "0" (m), "c" (n), "D" (f), "S" (g) : "%rbx"
	);
	return p;
}

uint64_t mpn_add_n (uint64_t *f, const uint64_t *g, const uint64_t *h, long n)
{
	uint64_t p;
	__asm__ __volatile__ (
		"xorq  %%rbx, %%rbx\n\t"
		"L1%=:\n\t"
		"movq (%%rsi, %%rbx, 8), %%rax\n\t"
		"adcq (%%rdx, %%rbx, 8), %%rax\n\t"
		"movq  %%rax,(%%rdi, %%rbx, 8)\n\t"
		"incq  %%rbx\n\t"
		"loop  L1%=\n\t"
		"movq  $1, %%rax\n\t"
		"jc    L2%=\n\t"
		"xorq  %%rax, %%rax\n\t"
		"L2%=:\n\t"
		: "=a" (p) : "c" (n), "d" (h), "S" (g), "D" (f) : "%rbx"
	);
	return p;
}

uint64_t mcn_invm(uint64_t a, uint64_t mod)
{
	uint64_t r;
	__asm__ __volatile__ (
		"movq  $0,    %%rdi # s = 0 \n\t"
		"movq  $0,    %%rsi # t = 0 \n\t"
		"movq  $0,    %%r8  # u = 0 \n\t"
		"movq  $1,    %%r9  # v = 1 \n\t"
		"movq  %%rcx, %%rax # mod -> rax \n\t"
	"L1%=:\n\t"
		"testq %%rbx, %%rbx #   a -> rbx \n\t"
		"jz  L6%=\n\t"
		"movq  $0,    %%rdx\n\t"
		"divq  %%rbx        # q = mod / a \n\t"
		"pushq %%rdx\n\t"
		"mulq  %%r9         # w = q * v \n\t"
		"popq  %%rdx\n\t"
		"cmpq  %%rax,  %%r8\n\t"
		"jnb L2%=\n\t"
		"xchgq %%r8,  %%rax\n\t"
		"movq  %%rdi, %%rsi\n\t"
		"notq  %%rdi\n\t"
	"L3%=:\n\t"
		"cmpq  %%rsi, %%rdi\n\t"
		"je  L4%=\n\t"
		"addq  %%r8,  %%rax\n\t"
	"L5%=:\n\t"
		"movq  %%r9,  %%r8  # v -> u \n\t"
		"movq  %%rax, %%r9  # u -> v \n\t"
		"movq  %%rbx, %%rax # a -> b \n\t"
		"movq  %%rdx, %%rbx # b -> a \n\t"
		"jmp L1%=\n\t"
	"L2%=:\n\t"
		"xchgq %%rdi, %%rsi\n\t"
		"jmp L3%=\n\t"
	"L4%=:\n\t"
		"subq  %%r8,  %%rax\n\t"
		"jmp L5%=\n\t"
	"L6%=:\n\t"
		"movq  %%r8, %%rax\n\t"
		"testq %%rsi, %%rsi\n\t"
		"jz  L7%=\n\t"
		"movq  %%rcx, %%rax\n\t"
		"subq  %%r8,  %%rax\n\t"
	"L7%=:\n\t"
		: "=a" (r) : "b"(a), "c"(mod) : "%rdi", "%rsi"
	);
	return r;
}

uint32_t mcn_invm(uint32_t a, uint32_t mod)
{
	uint32_t r;
	__asm__ __volatile__ (
		"movl  $0,    %%edi # s = 0 \n\t"
		"movl  $0,    %%esi # t = 0 \n\t"
		"movl  $0,    %%r8d # u = 0 \n\t"
		"movl  $1,    %%r9d # v = 1 \n\t"
		"movl  %%ecx, %%eax # mod -> eax \n\t"
	"L1%=:\n\t"
		"testl %%ebx, %%ebx #   a -> ebx \n\t"
		"jz  L6%=\n\t"
		"movl  $0,    %%edx\n\t"
		"divl  %%ebx        # q = mod / a \n\t"
		"pushq %%rdx\n\t"
		"mull  %%r9d        # w = q * v \n\t"
		"popq  %%rdx\n\t"
		"cmpl  %%eax, %%r8d\n\t"
		"jnb L2%=\n\t"
		"xchgl %%r8d, %%eax\n\t"
		"movl  %%edi, %%esi\n\t"
		"notl  %%edi\n\t"
	"L3%=:\n\t"
		"cmpl  %%esi, %%edi\n\t"
		"je  L4%=\n\t"
		"addl  %%r8d, %%eax\n\t"
	"L5%=:\n\t"
		"movl  %%r9d, %%r8d # v -> u \n\t"
		"movl  %%eax, %%r9d # u -> v \n\t"
		"movl  %%ebx, %%eax # a -> b \n\t"
		"movl  %%edx, %%ebx # b -> a \n\t"
		"jmp L1%=\n\t"
	"L2%=:\n\t"
		"xchgl %%edi, %%esi\n\t"
		"jmp L3%=\n\t"
	"L4%=:\n\t"
		"subl  %%r8d, %%eax\n\t"
		"jmp L5%=\n\t"
	"L6%=:\n\t"
		"movl  %%r8d, %%eax\n\t"
		"testl %%esi, %%esi\n\t"
		"jz  L7%=\n\t"
		"movl  %%ecx, %%eax\n\t"
		"subl  %%r8d, %%eax\n\t"
	"L7%=:\n\t"
		: "=a" (r) : "b"(a), "c"(mod) : "%edi", "%esi"
	);
	return r;
}
