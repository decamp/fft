package bits.fft;

import java.util.Map;
import java.util.WeakHashMap;


/**
 * This is a very simple convenience utility for timing events.
 * Just call <code>start()</code>, and then anytime you call
 * <code>seconds()</code> it will tell you have many seconds since
 * the calling thread called <code>start()</code>.
 * <p/>
 * Data is stored by calling Thread.  Different threads will
 * not interfere, although you may get unexpected results if
 * you try to time events that occur across multiple threads.
 *
 * @author Philip DeCamp
 */
public class Timer {

    private static Map<Thread, Long> mTimes = new WeakHashMap<Thread, Long>();


    /**
     * Start the timer.
     */
    public static synchronized void start() {
        Thread t = Thread.currentThread();
        mTimes.put( t, System.nanoTime() );
    }

    /**
     * Returns seconds since <code>start()</code> was last called by calling Thread,
     * or -1.0 if not called.
     *
     * @return seconds
     */
    public static synchronized double seconds() {
        long now = System.nanoTime();
        Long start = mTimes.remove( Thread.currentThread() );

        if( start == null ) {
            return -1.0;
        }

        return (now - start) / 1000000000.0;
    }

    /**
     * Same as <code>seconds()</code>, but prints message to System.out.
     *
     * @param msg
     * @return seconds
     */
    public static synchronized double printSeconds( String msg ) {
        double v = seconds();
        System.out.println( msg + v );
        return v;
    }

}
