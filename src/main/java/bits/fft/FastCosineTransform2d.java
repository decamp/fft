/*
 * Copyright (c) 2014, Massachusetts Institute of Technology
 * Released under the BSD 2-Clause License
 * http://opensource.org/licenses/BSD-2-Clause
 */
package bits.fft;

/**
 * Performs a Fast Discrete Cosine Transform on an square matrix of real values.
 * Currently relies on FastFourierTransform2d, so it could by faster
 * if it used a more specialized FCT algorithm.
 * <p/>
 * Not thread safe. For parallel operation, create multiple FastCosineTransform2d
 * objects.
 *
 * @author Philip DeCamp
 */
public class FastCosineTransform2d {


    /**
     * Factory method to create FastCosineTransform2d objects.
     * The sole argument, <code>dim</code>, indicates the size
     * of the matrices on which this transform will operate,
     * which must be a power-of-two.  For example, if dim == 8,
     * then all input arrays must be have a length of at least
     * 8*8.
     * <p/>
     * Memory allocation is just above (dim*dim*4+dim*4)*8 bytes.
     * <p/>
     *
     * @param dim Size of one-side of square matrix on which this transform operates.  Must be power-of-two.
     * @return a new FastCosineTransform2d instance.
     * @throws IllegalArgumentException if dim is smaller than 2, too large, or not a power-of-two.
     */
    public static FastCosineTransform2d create( int dim ) {
        return new FastCosineTransform2d( dim );
    }

    /**
     * Performs a Fast Discrete Cosine Transform on a square matrix of real values.
     * <br/>
     * Samples must be stored in a tight-packed format. For a 2x2 matrix: <br/>
     * <code>[... r_0_0, r_1_0, r_0_1, r_1_1, ...] </code><br/>,
     * where <code>r_m_n</code> is the real-valued element at position [m,n].
     *
     * @param a       Input matrix of real values with size [dim,dim].
     * @param aOff    Offset into array <code>a</code>
     * @param inverse Set to <code>true</code> to perform inverse transform.
     * @param out     Output matrix where DCT coeffs are stored.  Must have space for <code>dim*dim</code> values.
     * @param outOff  Offset into array<code>out</code>
     */
    public void apply( double[] a, int aOff, boolean inverse, double[] out, int outOff ) {
        if( !inverse ) {
            shuffle1( a, aOff, mDim, mBits, mWorkA );
            FastFourierTransform2d.doFft2d( mWorkA, 0, mDim, false );
            shuffle2( mWorkA, mWeight, mDim, mBits, mWorkB );
            FastFourierTransform2d.doFft2d( mWorkB, 0, mDim, false );
            shuffle3( mWorkB, mWeight, mDim, out, outOff );
        } else {
            invShuffle1( a, aOff, mInvWeight, mDim, mBits, mWorkA );
            FastFourierTransform2d.doFft2d( mWorkA, 0, mDim, true );
            invShuffle2( mWorkA, mInvWeight, mDim, mBits, mWorkB );
            FastFourierTransform2d.doFft2d( mWorkB, 0, mDim, true );
            invShuffle3( mWorkB, mDim, mBits, out, outOff );
        }
    }



    private static final int MAX_BITS = 15;

    private final int mDim;
    private final int mBits;

    private final double[] mWeight;
    private final double[] mInvWeight;
    private final double[] mWorkA;
    private final double[] mWorkB;


    private FastCosineTransform2d( int dim ) {
        mDim  = dim;
        mBits = computeNumBits( dim );

        mWeight    = new double[dim * 2];
        mInvWeight = new double[dim * 2];
        mWorkA     = new double[dim * dim * 2];
        mWorkB     = new double[dim * dim * 2];

        computeWeightVector( dim, mWeight );
        invComputeWeightVector( dim, mInvWeight );
    }


    private static int computeNumBits( int n ) {
        int bits = 1; //Start at 1 because 2^0 is invalid.
        for(; bits < MAX_BITS; bits++ ) {
            if( (1 << bits) == n ) {
                break;
            }
        }

        if( bits == MAX_BITS ) {
            throw new IllegalArgumentException( "Dimension must be a power of two, larger than 1, and smaller than 1 << " + MAX_BITS );
        }

        return bits;
    }


    private static void computeWeightVector( int dim, double[] w) {
        w[0] = 1.0;
        w[1] = 0.0;
        for( int i = 1; i < dim; i++ ) {

            double cos = Math.cos( -i * Math.PI * 0.5 / dim );
            double sin = -Math.sqrt( 1.0 - cos * cos );
            w[i*2  ] = 2.0 * cos;
            w[i*2+1] = 2.0 * sin;
        }
    }


    private static void invComputeWeightVector( int dim, double[] w ) {
        final double s = 1.0 / dim;
        w[0] = 1.0 * s;
        w[1] = 0.0;

        for( int i = 1; i < dim; i++ ) {
            double cos = Math.cos( i * Math.PI * 0.5 / dim );
            double sin = Math.sqrt( 1.0 - cos * cos );
            w[i*2  ] = cos * s;
            w[i*2+1] = sin * s;
        }
    }


    /**
     * 1. Drop complex components <br/>
     * 2. Butterfly shuffle <br/>
     * 3. Reverse-bit shuffle <br/>
     */
    private static void shuffle1( double[] a, int offA, int dim, int bits, double[] out ) {
        final int shift = 30 - bits;
        final int dim2 = dim * 2;

        for( int rowOut2 = 0; rowOut2 < dim2; rowOut2 += 2 ) {
            //Bit reversal shuffle.
            int rowIn = (reverseBits( rowOut2 ) >>> shift);

            //Butterfly shuffle.  
            //Equivalent to: if(x < dim / 2) {ia = x * 2;} else {ia = dim*2 - x*2 - 1;}
            //Example where dim = 8: [0 1 2 3 4 5 6 7] -> [0 2 4 6 7 5 3 1] 
            rowIn = rowIn + (rowIn / dim) * (dim2 - 2 * rowIn - 1);

            //Translate row A position into array element.
            int indA = rowIn + offA;

            for( int y = 0; y < dim; y++ ) {
                int vec = y * dim2;
                out[rowOut2 + vec] = a[indA + y * dim];
                out[rowOut2 + vec + 1] = 0.0;
            }
        }
    }

    /**
     * 1. Apply weight vector. <br/>
     * 2. Drow imaginary components. <br/>
     * 3. Transpose. <br/>
     * 4. Butterfly shuffle. <br/>
     * 5. Reverse-bit shuffle. <br/>.
     */
    private static void shuffle2( double[] a, double[] w, int dim, int bits, double[] out ) {
        final int shift = 30 - bits;
        final int dim2 = dim * 2;

        for( int rowOut2 = 0; rowOut2 < dim2; rowOut2 += 2 ) {
            //Bit reversal shuffle.
            int rowA = (reverseBits( rowOut2 ) >>> shift);

            //Butterfly shuffle.  
            //Equivalent to: if(x < dim / 2) {ia = x * 2;} else {ia = dim*2 - x*2 - 1;}
            //Example where dim = 8: [0 1 2 3 4 5 6 7] -> [0 2 4 6 7 5 3 1] 
            rowA = rowA + (rowA / dim) * (dim2 - 2 * rowA - 1);

            //Translate row A position into array offsets.
            //Note that A is being transposed.
            int indA = rowA * dim2;

            //Iterater through each column.
            for( int col2 = 0; col2 < dim2; col2 += 2 ) {
                //Perform complex multiplaction; throw away imaginary component.
                out[rowOut2 + col2 * dim] = a[indA + col2] * w[col2] - a[indA + col2 + 1] * w[col2 + 1];
                out[rowOut2 + col2 * dim + 1] = 0.0;
            }
        }
    }

    /**
     * 1. Apply weight vector.
     * 2. Transpose
     * 3. Drop imaginary components.
     *
     * @param a      Input: complex matrix of size [dim,dim]
     * @param w      Weight vector of size [dim] needed to compute DCT from FFT.
     * @param dim    Dimension of transform.  Must be power-of-two.
     * @param out    Output: real-valued matrix of size [dim,dim].
     * @param offOut Offset into <code>out</code>.
     */
    private static void shuffle3( double[] a, double[] w, int dim, double[] out, int offOut ) {
        final int dim2 = 2 * dim;

        for( int y2 = 0; y2 < dim2; y2 += 2 ) {
            double vReal = w[y2];
            double vImag = w[y2 + 1];

            for( int x = 0; x < dim; x++ ) {
                int ii = y2 + x * dim2;
                out[offOut++] = a[ii] * vReal - a[ii + 1] * vImag;
            }
        }
    }

    /**
     * 1. apply weights, <br/>
     * 2. reverse-bit shuffle
     */
    private static void invShuffle1( double[] a, int offA, double[] w, int dim, int bits, double[] out ) {
        final int dim2 = dim * 2;
        final int shift = 31 - bits;

        for( int rowIn = 0; rowIn < dim; rowIn++ ) {
            int indA = offA + rowIn;
            int indOut = (reverseBits( rowIn ) >>> shift);
            double wReal = w[rowIn * 2];
            double wImag = w[rowIn * 2 + 1];

            for( int col = 0; col < dim; col++ ) {
                double v = a[indA + col * dim];
                out[indOut + col * dim2] = wReal * v;
                out[indOut + col * dim2 + 1] = wImag * v;
            }
        }
    }

    /**
     * 1. drop imaginary component<br/>
     * 2. reverse butterfly shuffle<br/>
     * 3. transpose
     * 4. apply weights
     * 5. bit reversal shuffle<br/>
     */
    private static void invShuffle2( double[] a, double[] w, int dim, int bits, double[] out ) {
        final int shift = 30 - bits;
        final int dim2 = dim * 2;

        for( int i2 = 0; i2 < dim2; i2 += 2 ) {
            //Bit reverse shuffle.
            final int col2 = reverseBits( i2 ) >>> shift;
            final int indA = col2 * dim;

            final double wReal = w[col2];
            final double wImag = w[col2 + 1];

            for( int j2 = 0; j2 < dim2; j2 += 2 ) {
                //Reverse butterfly shuffle.
                int rowO = j2 + (j2 / dim) * (dim2 - 1 - 2 * j2);
                double val = a[j2 + indA];

                out[i2 + rowO * dim2] = val * wReal;
                out[i2 + rowO * dim2 + 1] = val * wImag;
            }
        }
    }

    /**
     * 1. drop imaginary component<br/>
     * 2. reverse butterfly shuffle<br/>
     */
    private static void invShuffle3( double[] a, int dim, int bits, double[] out, int offOut ) {
        final int shift = 32 - bits;
        final int dim2 = dim * 2;

        for( int i2 = 0; i2 < dim2; i2 += 2 ) {
            //Reverse butterfly shuffle.
            final int rowO = i2 + (i2 / dim) * (dim2 - 1 - 2 * i2);
            final int indO = offOut + rowO * dim;

            for( int j = 0; j < dim; j++ ) {
                out[j + indO] = a[i2 + j * dim2];
            }
        }
    }

    /**
     * @return Version of val with reversed bits.
     */
    private static int reverseBits( int val ) {
        long v = (((((val >>> 24)       ) * 0x0202020202L & 0x010884422010L) % 1023L)      ) |
                 (((((val >>> 16) & 0xFF) * 0x0202020202L & 0x010884422010L) % 1023L) << 8 ) |
                 (((((val >>>  8) & 0xFF) * 0x0202020202L & 0x010884422010L) % 1023L) << 16) |
                 (((((val       ) & 0xFF) * 0x0202020202L & 0x010884422010L) % 1023L) << 24);

        return (int)v;
    }

}