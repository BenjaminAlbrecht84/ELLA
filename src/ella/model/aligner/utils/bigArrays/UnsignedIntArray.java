package ella.model.aligner.utils.bigArrays;

import java.util.Arrays;

public class UnsignedIntArray {

    // TODO: apply http://www.omsn.de/blog/big-arrays-in-java

    private final static int MAX_SIZE = Integer.MAX_VALUE - 8;
    private int[] a1, a2;
    private int size;

    public UnsignedIntArray(UnsignedIntArray a) {
        this.size = (int) a.size();
        int[] a1Clone = (int[]) a.getArrays()[0];
        a1 = Arrays.copyOf(a1Clone, a1Clone.length);
        int[] a2Clone = (int[]) a.getArrays()[1];
        a2 = Arrays.copyOf(a2Clone, a2Clone.length);
    }

    public UnsignedIntArray(int[] a1, int[] a2) {
        this.size = a1.length + a2.length;
        this.a1 = a1;
        this.a2 = a2;
    }

    public UnsignedIntArray(UnsignedIntArray a, long size) {
        this.size = (int) size;
        int[] split = splitLong(size);
        if (size > a.size()) {
            int[] a1Clone = (int[]) a.getArrays()[0];
            a1 = new int[split[0]];
            System.arraycopy(a1Clone, 0, a1, 0, a1Clone.length);
            int[] a2Clone = (int[]) a.getArrays()[1];
            a2 = new int[split[1]];
            System.arraycopy(a2Clone, 0, a2, 0, a2Clone.length);
        } else {
            int[] a1Clone = (int[]) a.getArrays()[0];
            a1 = new int[split[0]];
            System.arraycopy(a1Clone, 0, a1, 0, (int) (size < MAX_SIZE ? size : MAX_SIZE));
            int[] a2Clone = (int[]) a.getArrays()[1];
            a2 = new int[split[1]];
            System.arraycopy(a2Clone, 0, a2, 0, (int) (size > MAX_SIZE ? size - MAX_SIZE : 0));
        }
    }

    public UnsignedIntArray(long size) {
        this.size = (int) size;
        int[] split = splitLong(size);
        a1 = new int[split[0]];
        a2 = new int[split[1]];
    }

    public long get(long index) {
        int[] split = splitLong(index);
        int value = split[0] < a1.length ? a1[split[0]] :  a2[split[1]];
        return Integer.toUnsignedLong(value);
    }

    public void set(long index, long value) {
        int[] split = splitLong(index);
        if (split[0] < a1.length)
            a1[split[0]] = (int) value;
        else if ((split[1] < a2.length))
            a2[split[1]] = (int) value;
    }

    public long size() {
        return Integer.toUnsignedLong(size);
    }

    private int[] splitLong(long value) {
        int v1, v2;
        if (value < MAX_SIZE) {
            v1 = (int) value;
            v2 = 0;
        } else {
            v1 = MAX_SIZE;
            v2 = (int) (value - MAX_SIZE);
        }
        int[] result = {v1, v2};
        return result;
    }

    public Object[] getArrays() {
        Object[] arrays = {a1, a2};
        return arrays;
    }

}
