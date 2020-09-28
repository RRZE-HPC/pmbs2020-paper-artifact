#define INSTR ld1d
#define NINST 16
#define STRIDE #0
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

        mov     z16.d, p0/z, STRIDE
        mov     z17.d, p0/z, STRIDE
        mov     z18.d, p0/z, STRIDE
        mov     z19.d, p0/z, STRIDE
        mov     z20.d, p0/z, STRIDE
        mov     z21.d, p0/z, STRIDE
        mov     z22.d, p0/z, STRIDE
        mov     z23.d, p0/z, STRIDE
        mov     z24.d, p0/z, STRIDE
        mov     z25.d, p0/z, STRIDE
        mov     z26.d, p0/z, STRIDE
        mov     z27.d, p0/z, STRIDE
        mov     z28.d, p0/z, STRIDE
        mov     z29.d, p0/z, STRIDE
        mov     z30.d, p0/z, STRIDE
        mov     z31.d, p0/z, STRIDE
    
        mov     x1, sp
        # align to 64 byte
        and     x1, x1, #0xFFFFFFFFFFFFFFC0

loop:
        subs      x4, x4, #1
        INSTR    z0.d, p0/z, [x1, z20.d, lsl 3]
        INSTR    z1.d, p0/z, [x1, z21.d, lsl 3]
        INSTR    z2.d, p0/z, [x1, z22.d, lsl 3]
        INSTR    z3.d, p0/z, [x1, z23.d, lsl 3]
        INSTR    z4.d, p0/z, [x1, z24.d, lsl 3]
        INSTR    z5.d, p0/z, [x1, z25.d, lsl 3] 
        INSTR    z6.d, p0/z, [x1, z26.d, lsl 3]
        INSTR    z7.d, p0/z, [x1, z27.d, lsl 3] 
        INSTR    z8.d, p0/z, [x1, z28.d, lsl 3]
        INSTR    z9.d, p0/z, [x1, z29.d, lsl 3] 
        INSTR    z10.d, p0/z, [x1, z30.d, lsl 3]
        INSTR    z11.d, p0/z, [x1, z31.d, lsl 3] 
        INSTR    z12.d, p0/z, [x1, z16.d, lsl 3]
        INSTR    z13.d, p0/z, [x1, z17.d, lsl 3] 
        INSTR    z14.d, p0/z, [x1, z18.d, lsl 3]
        INSTR    z15.d, p0/z, [x1, z19.d, lsl 3] 
        bne       loop
done:
        # pop callee-save registers from stack
        ld1            {v12.2d, v13.2d, v14.2d, v15.2d}, [sp]
        add            sp, sp, #64
        ld1            {v8.2d, v9.2d, v10.2d, v11.2d}, [sp]
        add            sp, sp, #64
        
        ret

.size latency, .-latency
