#define INSTR ld1d
#define NINST 16
#define STRIDE1 #6
#define STRIDE2 #8
#define N x0

.arch armv8.2-a+sve
.globl ninst
.data
ninst:
.long NINST
.text
.globl latency
.type latency, @function
.align 2
latency:

        # push callee-save registers onto stack
        sub            sp, sp, #64
        st1            {v8.2d, v9.2d, v10.2d, v11.2d}, [sp]
        sub            sp, sp, #64
        st1            {v12.2d, v13.2d, v14.2d, v15.2d}, [sp]

        mov     x4, N
        
        ptrue   p0.d
        ptrue   p2.d, vl8
        pfalse  p1.b
        trn1    p1.d, p2.d, p1.d

        mov     z0.d, p0/z, STRIDE1
        mov     z0.d, p1/m, STRIDE2

        mov     x1, sp
        # align to 64 byte
        and     x1, x1, #0xFFFFFFFFFFFFFFC0

loop:
        subs      x4, x4, #1
        INSTR    z12.d, p0/z, [x1, z0.d, lsl 3]
        INSTR    z13.d, p0/z, [x1, z0.d, lsl 3] 
        INSTR    z14.d, p0/z, [x1, z0.d, lsl 3]
        INSTR    z15.d, p0/z, [x1, z0.d, lsl 3] 
        INSTR    z16.d, p0/z, [x1, z0.d, lsl 3]
        INSTR    z1.d,  p0/z, [x1, z0.d, lsl 3]
        INSTR    z2.d,  p0/z, [x1, z0.d, lsl 3]
        INSTR    z3.d,  p0/z, [x1, z0.d, lsl 3]
        INSTR    z4.d,  p0/z, [x1, z0.d, lsl 3]
        INSTR    z5.d,  p0/z, [x1, z0.d, lsl 3] 
        INSTR    z6.d,  p0/z, [x1, z0.d, lsl 3]
        INSTR    z7.d,  p0/z, [x1, z0.d, lsl 3] 
        INSTR    z8.d,  p0/z, [x1, z0.d, lsl 3]
        INSTR    z9.d,  p0/z, [x1, z0.d, lsl 3] 
        INSTR    z10.d, p0/z, [x1, z0.d, lsl 3]
        INSTR    z11.d, p0/z, [x1, z0.d, lsl 3] 
        bne       loop
done:
        # pop callee-save registers from stack
        ld1            {v12.2d, v13.2d, v14.2d, v15.2d}, [sp]
        add            sp, sp, #64
        ld1            {v8.2d, v9.2d, v10.2d, v11.2d}, [sp]
        add            sp, sp, #64
        
        ret

.size latency, .-latency
