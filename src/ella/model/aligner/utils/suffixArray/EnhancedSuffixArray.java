package ella.model.aligner.utils.suffixArray;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;

import ella.model.aligner.aligning.Aligner;
import ella.model.aligner.aligning.query.QueryContainer;
import ella.model.aligner.indexing.IndexText;
import ella.model.aligner.utils.Alphabet;
import ella.model.aligner.utils.CyclicSeedShape;
import ella.model.aligner.utils.bigArrays.UnsignedIntArray;
import ella.model.aligner.utils.streams.MyIntegerStream;
import ella.model.aligner.utils.wrapper.ByteArrayWrapper;

public class EnhancedSuffixArray {

    private ByteArrayWrapper commonPrefix;
    private CyclicSeedShape seedShape;

    private UnsignedIntArray suffixArray;
    private LCPTreeArray lcpTreeArray;
    private BucketTable bucketTable;
    private SearchTable searchTable;

    public EnhancedSuffixArray(ByteArrayWrapper prefix, CyclicSeedShape seedShape) {
        this.commonPrefix = prefix;
        this.seedShape = seedShape;
    }

    public EnhancedSuffixArray(ByteArrayWrapper prefix, CyclicSeedShape seedShape, UnsignedIntArray suffixArray, LCPTreeArray lcpTreeArray,
                               BucketTable bucketTable) {
        this.commonPrefix = prefix;
        this.seedShape = seedShape;
        this.suffixArray = suffixArray;
        this.lcpTreeArray = lcpTreeArray;
        this.bucketTable = bucketTable;
        // this.searchTable = new SearchTable(seedShape, prefix.getAAPrefix());
    }

    public void initSuffixArray(IndexText text, UnsignedIntArray suffixes, MSDRadixSortForSuffixes radixSort) {
        Object[] result = radixSort.run(text, suffixes, seedShape, commonPrefix);
        suffixArray = (UnsignedIntArray) result[0];
        lcpTreeArray = (LCPTreeArray) result[1];
        bucketTable = (BucketTable) result[2];
    }

    // computing the rightmost prefix having less than or equal to 'm' hits
    public ArrayList<Object[]> findPattern(IndexText text, String[] splitPattern, int minLen, int m, boolean useBucketTable, QueryContainer qC,
                                           Aligner aligner, ImprovedBinarySearch search) {

        ArrayList<Object[]> matches;
        if (suffixArray.size() <= m) {
            matches = new ArrayList<>((int) suffixArray.size());
            for (int i = 0; i < suffixArray.size(); i++) {
                Object[] match = {suffixArray.get(i), commonPrefix.getLength(), i};
                matches.add(match);
            }
            if (!splitPattern[0].isEmpty())
                checkMatches(text, matches, splitPattern[1]);
        } else {
            // int[] range = searchTable.findBucketEntry(splitPattern[1]);
            // if (range != null)
            // return cmpRangeResult(range, splitPattern);
//            matches = new ImprovedBinarySearch_2().run(text, splitPattern, this, seedShape, minLen, m, null, commonPrefix, qC, aligner);
            matches = search.run(text, splitPattern, this, seedShape, minLen, m, null, commonPrefix, qC, aligner);

        }

        return matches;
    }

    private ArrayList<Object[]> cmpRangeResult(int[] range, String[] splitPattern) {
        ArrayList<Object[]> matches = new ArrayList<Object[]>();
        long l = Integer.toUnsignedLong(range[0]);
        long r = Integer.toUnsignedLong(range[1]);
        for (long i = l; i <= r; i++) {
            Object[] match = {suffixArray.get(i), range[2] + splitPattern[0].length(), i};
            matches.add(match);
        }
        return matches;
    }

    public void updateSearchTable(String s, long[] range, int length) {
        // searchTable.setBucketEntry(s, range, length);
    }

    private void checkMatches(IndexText text, ArrayList<Object[]> matches, String header) {
        MyIntegerStream deleteStream = new MyIntegerStream(matches.size());
        for (int n = 0; n < matches.size(); n++) {
            long index = (long) matches.get(n)[0];
            int pos = 0;
            for (int i = 1; i <= header.length(); i++) {

                // protein beginning reached, no match
                if ((index - i) < 0 || text.readPosition(index - i) == 127) {
                    deleteStream.add(n);
                    break;
                }

                // checking characters
                if (!seedShape.usePosition(pos++))
                    continue;
                char c1 = Alphabet.reduceCharacter(header.charAt(header.length() - i));
                char c2 = Alphabet.reduceCharacter(text.readAA(index - i));
                if (c1 != c2) {
                    deleteStream.add(n);
                    break;
                }

            }
        }
        UnsignedIntArray toDelete = deleteStream.toArray();
        for (long i = toDelete.size() - 1; i >= 0; i--)
            matches.remove(toDelete.get(i));

    }

    public ByteArrayWrapper getCommonPrefix() {
        return commonPrefix;
    }

    public CyclicSeedShape getSeedShape() {
        return seedShape;
    }

    public UnsignedIntArray getSuffixArray() {
        return suffixArray;
    }

    public LCPTreeArray getLcpTreeArray() {
        return lcpTreeArray;
    }

    public BucketTable getBucketTable() {
        return bucketTable;
    }

    public String getAAPrefix() {
        StringBuilder builder = new StringBuilder();
        for (byte b : commonPrefix.getData()) {
            byte[] bytes = {b};
            ByteBuffer buffer = ByteBuffer.wrap(bytes);
            buffer.order(ByteOrder.LITTLE_ENDIAN);
            builder.append(Alphabet.getCharacter((int) buffer.get()));
        }
        return builder.toString();
    }

    // ONLY FOR DEBUGGING

    public void printSuffixArray(IndexText text) {
        if (suffixArray != null) {
            System.out.println("---");
            System.out.println(getAAPrefix() + " " + seedShape.toString());
            UnsignedIntArray array = suffixArray;
            for (long i = 0; i < suffixArray.size(); i++) {
                int lcp = lcpTreeArray == null ? -1 : lcpTreeArray.getTreeValue(i - 1, i);
                System.out.print(i + "\t" + array.get(i) + " \t " + lcp + "\t");
                for (int j = 0; j < 50; j++) {
                    int index = text.readPosition(array.get(i) + j);
                    if (index == 127)
                        break;
                    System.out.print(Alphabet.getReducedCharacter(index));
                }
                System.out.println();
            }
        }
    }

    public long getByteSize() {
        if (suffixArray == null)
            return 0;
        return suffixArray.size() * 4 + lcpTreeArray.getIntSize() * 4;
    }

    public long getSAByteSize() {
        if (suffixArray == null)
            return 0;
        return suffixArray.size() * 4;
    }

    public long getLCPByteSize() {
        if (lcpTreeArray == null)
            return 0;
        return lcpTreeArray.getIntSize() * 4;
    }

    public void freeMemory() {
        suffixArray = null;
        lcpTreeArray.freeMemory();
    }

}
