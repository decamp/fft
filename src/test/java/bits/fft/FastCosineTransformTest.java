/*
 * Copyright (c) 2014, Massachusetts Institute of Technology
 * Released under the BSD 2-Clause License
 * http://opensource.org/licenses/BSD-2-Clause
 */
package bits.fft;

import org.junit.Test;
import java.util.Random;

public class FastCosineTransformTest {


    private static final int OFF_0 = 3;
    private static final int DIM_0 = 16;

    private static final double[] INPUT_0 = {
            0.0,
            0.0,
            0.0,
            0.4617563814065817,
            -0.17983837701559668,
            -0.5845703173805659,
            -0.33456588808097765,
            0.9355118188482414,
            -0.9877656354684774,
            0.9274095940464153,
            0.8797307775638197,
            0.8943898353263877,
            0.8741642977919393,
            -0.2056513156305888,
            -0.3049639415937795,
            -0.41188593599192647,
            0.012967254652470173,
            -0.7680658239346845,
            0.5410717601583555
    };

    private static final double[] OUTPUT_0 = {
            1.749694484697614,
            0.8980720787935931,
            -5.931515629795895,
            -0.30061496357145784,
            7.349553131552409,
            -0.8691238683175415,
            0.1004889009271066,
            2.4211891107914396,
            5.052129643815888,
            0.6438664036160185,
            2.489033936181609,
            -5.505819456960973,
            -2.6246140382688674,
            -1.3605822994598682,
            6.48923375792617,
            2.409897519888083
    };


    @Test
    public void testCorrect() {
        final int dim = DIM_0;
        final int off = OFF_0;

        double[] b = new double[dim + off];
        double[] c = new double[dim + off];

        FastCosineTransform trans = FastCosineTransform.create( dim );
        trans.apply( INPUT_0, off, false, b, off );
        trans.apply( b, off, true, c, off );

        //System.out.println(TestUtil.realToMatlab(INPUT_0, off, dim, 1));
        //System.out.println("===");
        //System.out.println(TestUtil.formatReal(b, off, dim, 1));
        //System.out.println(TestUtil.arrayToJava(b, off, dim, "OUTPUT_0"));


        TestUtil.assertNear( OUTPUT_0, 0, b, off, dim );
        TestUtil.assertNear( INPUT_0, off, c, off, dim );
    }


    @Test
    public void testSpeed() {
        final int dim = 512;

        double[] a = new double[dim];
        Random rand = new Random( 0 );

        for( int i = 0; i < dim; i++ ) {
            a[i] = rand.nextDouble() * 2.0 - 1.0;
        }


        double[] b = new double[dim];
        FastCosineTransform trans = FastCosineTransform.create( dim );

        Timer.start();

        for( int i = 0; i < 5000; i++ ) {
            trans.apply( a, 0, false, b, 0 );
        }

        Timer.printSeconds( "CosineTransform seconds: " );
    }

}