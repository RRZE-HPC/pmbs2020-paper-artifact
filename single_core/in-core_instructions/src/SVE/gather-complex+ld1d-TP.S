#define INSTR ld1d
#define NINST 10
#define STRIDE1 #0
#define STRIDE2 #16
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

        fmov    v0.2d, #1.00000000
        fmov    v1.2d, #1.00000000
        fmov    v2.2d, #1.00000000

        ptrue   p0.d
        ptrue   p2.d, vl8
        pfalse  p1.b
        trn1    p1.d, p2.d, p1.d

        mov     z20.d, p0/z, STRIDE1
        mov     z20.d, p1/m, STRIDE2
    
        mov     x0, sp
        # align to 64 byte
        and     x0, x0, #0xFFFFFFFFFFFFFFC0
        sub     x6, x1, #8192
        
        ldr     x1, [sp]
        st1d   {z20.d}, p0, [sp]
loop:
        subs      x4, x4, #1
        ld1d    z20.d, p1/z, [sp]
        INSTR   z0.d, p0/z, [x0, z20.d, lsl 3]
        ld1d    z21.d, p1/z, [sp]
        INSTR   z1.d, p0/z, [x0, z21.d, lsl 3]
        ld1d    z22.d, p1/z, [sp]
        INSTR   z2.d, p0/z, [x0, z22.d, lsl 3]
        ld1d    z23.d, p1/z, [sp]
        INSTR   z3.d, p0/z, [x0, z23.d, lsl 3]
        ld1d    z24.d, p1/z, [sp]
        INSTR   z4.d, p0/z, [x0, z24.d, lsl 3]
        ld1d    z25.d, p1/z, [sp]
        INSTR   z5.d, p0/z, [x0, z25.d, lsl 3] 
        ld1d    z26.d, p1/z, [sp]
        INSTR   z6.d, p0/z, [x0, z26.d, lsl 3]
        ld1d    z27.d, p1/z, [sp]
        INSTR   z7.d, p0/z, [x0, z27.d, lsl 3] 
        ld1d    z28.d, p1/z, [sp]
        INSTR   z8.d, p0/z, [x0, z28.d, lsl 3]
        ld1d    z29.d, p1/z, [sp]
        INSTR   z9.d, p0/z, [x0, z29.d, lsl 3] 
        bne       loop
done:
        str     x1, [sp]
        # pop callee-save registers from stack
        ld1            {v12.2d, v13.2d, v14.2d, v15.2d}, [sp]
        add            sp, sp, #64
        ld1            {v8.2d, v9.2d, v10.2d, v11.2d}, [sp]
        add            sp, sp, #64
        
        ret

.size latency, .-latency
