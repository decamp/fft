/*
 * Copyright (c) 2014, Massachusetts Institute of Technology
 * Released under the BSD 2-Clause License
 * http://opensource.org/licenses/BSD-2-Clause
 */
package bits.fft;

import static org.junit.Assert.assertTrue;


/**
 * @author Philip DeCamp
 */
class TestUtil {

    static String formatComplex( double[] v, int off, int rows, int cols ) {
        StringBuilder s = new StringBuilder();

        for( int y = 0; y < rows; y++ ) {
            for( int x = 0; x < cols; x++ ) {
                if( x > 0 ) {
                    s.append( "  " );
                }

                s.append( String.format( "% .4f %+.4f", v[off + (y + x * rows) * 2], v[off + (y + x * rows) * 2 + 1] ) );
            }
            s.append( '\n' );
        }

        return s.toString();
    }


    static String formatReal( double[] v, int off, int rows, int cols ) {
        StringBuilder s = new StringBuilder();

        for( int y = 0; y < rows; y++ ) {
            for( int x = 0; x < cols; x++ ) {
                if( x > 0 ) {
                    s.append( "  " );
                }

                s.append( String.format( "% .4f", v[off + (y + x * rows)] ) );
            }
            s.append( '\n' );
        }

        return s.toString();
    }


    static String complexToMatlab( double[] v, int off, int rows, int cols ) {
        StringBuilder s = new StringBuilder();
        s.append( "data = [" );

        for( int y = 0; y < rows; y++ ) {
            if( y > 0 ) {
                s.append( ';' );
            }

            for( int x = 0; x < cols; x++ ) {
                if( x > 0 ) {
                    s.append( ',' );
                }

                double vr = v[off + (x * rows + y) * 2];
                double vi = v[off + (x * rows + y) * 2 + 1];

                s.append( vr );

                if( vi >= 0.0 ) {
                    s.append( '+' );
                }

                s.append( vi ).append( 'i' );
            }
        }

        s.append( "];" );
        return s.toString();
    }


    static String realToMatlab( double[] v, int off, int rows, int cols ) {
        StringBuilder s = new StringBuilder();
        s.append( "data = [" );

        for( int y = 0; y < rows; y++ ) {
            if( y > 0 ) {
                s.append( ';' );
            }

            for( int x = 0; x < cols; x++ ) {
                if( x > 0 ) {
                    s.append( ',' );
                }

                s.append( v[off + (x * rows + y)] );
            }
        }

        s.append( "];" );
        return s.toString();
    }


    static String arrayToJava( double[] v, int off, int len, String name ) {
        String dec = "    private static final double[] " + name + " = { ";

        StringBuilder s = new StringBuilder();
        for( int i = 0; i < dec.length(); i++ ) {
            s.append( ' ' );
        }

        String prefix = s.toString();
        s = new StringBuilder( dec );

        for( int i = 0; i < len; i++ ) {
            if( i > 0 ) {
                s.append( ",\n" ).append( prefix );
            }
            s.append( v[i + off] );
        }

        s.append( "\n    };" );

        return s.toString();
    }


    static void assertNear( double[] a, int offA, double[] b, int offB, int n ) {
        for( int i = 0; i < n; i++ ) {
            double err = Math.abs( b[offB + i] - a[offA + i] );
            assertTrue( err < 0.0000001 );
        }
    }

}
