package ella.model.aligner.utils.suffixArray;

import java.util.ArrayList;
import java.util.stream.LongStream;

import ella.model.aligner.aligning.Aligner;
import ella.model.aligner.aligning.query.QueryContainer;
import ella.model.aligner.indexing.IndexText;
import ella.model.aligner.utils.Alphabet;
import ella.model.aligner.utils.CyclicSeedShape;
import ella.model.aligner.utils.bigArrays.UnsignedIntArray;
import ella.model.aligner.utils.streams.MyLongStream;
import ella.model.aligner.utils.wrapper.ByteArrayWrapper;

public class ImprovedBinarySearch {

    private static String query = "Escherichia-coli-str.-K-12-substr.-MG1655,-synthetic-genome_4174425_aligned_840_R_0_12711_0";
    private static String reference = "NP_418414.1";

    private long lastL = -1, lastR = -1;

    public ArrayList<Object[]> run(IndexText text, String[] splitPattern, EnhancedSuffixArray esa, CyclicSeedShape seedShape, int minLen, int m,
                                   long[] searchInterval, ByteArrayWrapper commonPrefix, QueryContainer qC, Aligner aligner) {

        boolean debug = false;

        // initializing fields
        UnsignedIntArray suffixArray = esa.getSuffixArray();
        LCPTreeArray lcpTreeArray = esa.getLcpTreeArray();
        if (suffixArray == null || lcpTreeArray == null)
            return null;

        // locating the position of the longest match by an improved binary search with O(m + log n) worst case running time
        Object[] longestMatch = binarySearch(splitPattern[1], text, suffixArray, lcpTreeArray, seedShape, searchInterval, m, commonPrefix);

        // checking neighboring matches of the longest match
        ArrayList<Object[]> allMatches = detectAllMatches(longestMatch, suffixArray, lcpTreeArray, minLen, m, debug, commonPrefix, qC, text,
                splitPattern);

        // updating search table
        updateSearchTable(allMatches, (int) longestMatch[1], splitPattern[1], esa);

        return allMatches;

    }

    private void updateSearchTable(ArrayList<Object[]> allMatches, int longestMatch, String pattern, EnhancedSuffixArray esa) {
        if (!allMatches.isEmpty() && pattern.length() > longestMatch) {
            long l = (long) allMatches.get(0)[2];
            long r = (long) allMatches.get(allMatches.size() - 1)[2];
            long[] range = {l, r};
            esa.updateSearchTable(pattern, range, longestMatch + 1);
        }
    }

    private ArrayList<Object[]> detectAllMatches(Object[] longestMatch, UnsignedIntArray suffixArray, LCPTreeArray lcpTreeArray, int minLen, int m,
                                                 boolean debug, ByteArrayWrapper commonPrefix, QueryContainer qC, IndexText text, String[] splitPattern) {

        ArrayList<Object[]> matches = new ArrayList<>();
        if ((int) longestMatch[1] < minLen || (int) longestMatch[1] == commonPrefix.getLength())
            return matches;

        if (splitPattern[0].isEmpty() || checkMatch(text, suffixArray.get((long) longestMatch[0]), splitPattern[0])) {
            Object[] match = {suffixArray.get((long) longestMatch[0]), (int) longestMatch[1] + splitPattern[0].length(), longestMatch[0]};
            matches.add(match);
        }

        // collecting all exact matches beginning with the longest ones (stopping as soon as there are more than 'm' exact matches)
        long l = (long) longestMatch[0];
        long r = (long) longestMatch[0] + 1;
        for (int len = (int) longestMatch[1]; len >= minLen; len--) {

            // finding exact matches of length 'len' on the left
            Object[] leftResult = findMatchesOnTheLeft(len, l, r, matches.size(), m, lcpTreeArray, suffixArray, splitPattern[0], text);
            MyLongStream leftMatches = (MyLongStream) leftResult[0];
            l = (long) leftResult[1];

            if (matches.size() + leftMatches.size() > m)
                break;

            // finding exact matches of length 'len' on the right
            Object[] rightResult = findMatchesOnTheRight(len, l, r, matches.size() + leftMatches.size(), m, lcpTreeArray, suffixArray,
                    splitPattern[0], text);
            MyLongStream rightMatches = (MyLongStream) rightResult[0];
            r = (long) rightResult[1];

            if (matches.size() + leftMatches.size() + rightMatches.size() > m)
                break;

            // processing lcp information
            addMatches(leftMatches, rightMatches, longestMatch, len, suffixArray, matches, text, splitPattern);

            // only taking the longest exact matches into account (this has been proven to be the best strategy)
            if (!matches.isEmpty())
                break;

        }

        return matches;

    }

    private void addMatches(MyLongStream leftMatches, MyLongStream rightMatches, Object[] longestMatch, int len, UnsignedIntArray suffixArray,
                            ArrayList<Object[]> matches, IndexText text, String[] splitPattern) {
        for (int i = 0; i < leftMatches.size(); i++) {
            long l = leftMatches.get(i);
            Object[] match = {suffixArray.get(l), len + splitPattern[0].length(), l};
            matches.add(0, match);
        }
        int pos = matches.size();
        for (int i = rightMatches.size() - 1; i >= 0; i--) {
            long r = rightMatches.get(i);
            Object[] match = {suffixArray.get(r), len + splitPattern[0].length(), r};
            matches.add(pos, match);
        }
    }

    private final Object[] leftMatchesResult = new Object[2];
    private final MyLongStream leftMatches = new MyLongStream();

    private Object[] findMatchesOnTheLeft(int len, long l, long r, int found, int m, LCPTreeArray lcpTreeArray, UnsignedIntArray suffixArray,
                                          String prefix, IndexText text) {
//        MyLongStream leftMatches = new MyLongStream(m);
        leftMatches.reset(m);
        while (l > 0 && found + leftMatches.size() <= m && r - l + 1 <= 5 * m) {
            int lcp = len < 128 ? lcpTreeArray.getUpperBoundValueAtPosition(l) : lcpTreeArray.getValueAtPosition(l);
            int length = Math.min(lcp, len);
            if (length < len)
                break;
            if (prefix.isEmpty() || checkMatch(text, suffixArray.get(l - 1), prefix))
                leftMatches.add(l - 1);
            l--;
        }
//        Object[] result = {leftMatches, l};
//        return result;
        leftMatchesResult[0] = leftMatches;
        leftMatchesResult[1] = l;
        return leftMatchesResult;
    }

    private final Object[] rightMatchesResult = new Object[2];
    private final MyLongStream rightMatches = new MyLongStream();

    private Object[] findMatchesOnTheRight(int len, long l, long r, int found, int m, LCPTreeArray lcpTreeArray, UnsignedIntArray suffixArray,
                                           String prefix, IndexText text) {
//        MyLongStream rightMatches = new MyLongStream(m);
        rightMatches.reset(m);
        while (r < suffixArray.size() && found + rightMatches.size() <= m && r - l + 1 <= 2 * m) {
            int lcp = len < 128 ? lcpTreeArray.getUpperBoundValueAtPosition(r) : lcpTreeArray.getValueAtPosition(r);
            int length = Math.min(lcp, len);
            if (length < len) {
                break;
            }
            if (prefix.isEmpty() || checkMatch(text, suffixArray.get(r), prefix))
                rightMatches.add(r);
            r++;
        }
//        Object[] result = {rightMatches, r};
//        return result;
        rightMatchesResult[0] = rightMatches;
        rightMatchesResult[1] = r;
        return rightMatchesResult;
    }

    private boolean checkMatch(IndexText text, long index, String prefix) {
        for (int i = 1; i <= prefix.length(); i++) {

            // checking if beginning of protein has been reached, no match in that case
            int pos = text.readPosition(index - i);
            if ((index - i) < 0 || pos == 127)
                return false;

            // checking for equal characters
            char c1 = Alphabet.reduceCharacter(prefix.charAt(prefix.length() - i));
            char c2 = Alphabet.getReducedCharacter(pos);
            if (c1 != c2)
                return false;

        }
        return true;
    }

    private Object[] binarySearch(String p, IndexText text, UnsignedIntArray suffixArray, LCPTreeArray lcpArray, CyclicSeedShape seedShape,
                                  long[] searchInterval, int m, ByteArrayWrapper commonPrefix) {

        long L = 0, R = suffixArray.size() - 1;
        if (searchInterval != null) {
            L = searchInterval[0];
            R = searchInterval[1] - 1;
        }
        int[] compareResult = compare(commonPrefix.getLength(), p, suffixArray.get(R), text, seedShape);

        int l = 0, r = compareResult[0];
        if (compareResult[1] < 0) {
            return binarySearchRec(L, divide(L + R, 2), R, l, r, p, text, suffixArray, lcpArray, seedShape, m, searchInterval == null);
        } else {
            Object[] result = {R, r};
            return result;
        }

    }

    private Object[] binarySearchRec(long L, long M, long R, int l, int r, String p, IndexText text, UnsignedIntArray suffixArray,
                                     LCPTreeArray lcpArray, CyclicSeedShape seedShape, int m, boolean useLCPTree) {

        if (L == R)
            l = r = Math.max(l, r);

        // reporting longest match if the entire pattern is not present
        if (L == lastL && R == lastR)
            return returnResult(L, M, R, l, r);
        lastL = L;
        lastR = R;

        // searching for pattern in interval [L,R]
        if (!useLCPTree || l == r) {
            return doTextComparisons(l, L, M, R, l, r, p, text, suffixArray, lcpArray, seedShape, useLCPTree, m);
        } else if (l > r) {
            // based on the LCP of L and M, either turning right, turning left, or doing a text comparison
            return caseOne(L, M, R, l, r, p, text, suffixArray, lcpArray, seedShape, useLCPTree, m);
        } else {
            // based on the LCP of R and M, either turning left, turning right, or doing a text comparison
            return caseTwo(L, M, R, l, r, p, text, suffixArray, lcpArray, seedShape, useLCPTree, m);
        }

    }

    private final Object[] binarySearchResult = new Object[2];

    private Object[] returnResult(long L, long M, long R, int l, int r) {
        if (l > r) {
//            Object[] result = {L, l};
//            return result;
            binarySearchResult[0] = L;
            binarySearchResult[1] = l;
            return binarySearchResult;
        } else {
//            Object[] result = {R, r};
//            return result;
            binarySearchResult[0] = R;
            binarySearchResult[1] = r;
            return binarySearchResult;
        }
    }

    private Object[] caseOne(long L, long M, long R, int l, int r, String p, IndexText text, UnsignedIntArray suffixArray, LCPTreeArray lcpArray,
                             CyclicSeedShape seedShape, boolean useLCPTree, int m) {

        int lcpLM = lcpArray.getTreeValue(L, M);

        if (lcpLM > l) {
            return binarySearchRec(M, divide(R + M, 2), R, l, r, p, text, suffixArray, lcpArray, seedShape, m, useLCPTree);
        } else if (lcpLM < l) {
            return binarySearchRec(L, divide(L + M, 2), M, l, lcpLM, p, text, suffixArray, lcpArray, seedShape, m, useLCPTree);
        } else {
            return doTextComparisons(l, L, M, R, l, r, p, text, suffixArray, lcpArray, seedShape, useLCPTree, m);
        }

    }

    private Object[] caseTwo(long L, long M, long R, int l, int r, String p, IndexText text, UnsignedIntArray suffixArray, LCPTreeArray lcpArray,
                             CyclicSeedShape seedShape, boolean useLCPTree, int m) {

        int lcpRM = lcpArray.getTreeValue(M, R);

        if (lcpRM > r) {
            return binarySearchRec(L, divide(L + M, 2), M, l, r, p, text, suffixArray, lcpArray, seedShape, m, useLCPTree);
        } else if (lcpRM < r) {
            return binarySearchRec(M, divide(R + M, 2), R, lcpRM, r, p, text, suffixArray, lcpArray, seedShape, m, useLCPTree);
        } else {
            return doTextComparisons(r, L, M, R, l, r, p, text, suffixArray, lcpArray, seedShape, useLCPTree, m);
        }
    }

    private final Object[] textComparisonResult = new Object[2];

    private Object[] doTextComparisons(int m, long L, long M, long R, int l, int r, String p, IndexText text, UnsignedIntArray suffixArray,
                                       LCPTreeArray lcpArray, CyclicSeedShape seedShape, boolean useLCPTree, int multiplicity) {

        // comparing pattern against a suffix of the text
        m = useLCPTree ? m : 0;
        int[] compareResult = compare(m, p, suffixArray.get(M), text, seedShape);
        m = compareResult[0];
        int comp = compareResult[1];

        // either reporting a match, turning left, or turning right
        if (comp == 0) {
//            Object[] result = {M, p.length()};
//            return result;
            textComparisonResult[0] = M;
            textComparisonResult[1] = p.length();
            return textComparisonResult;
        } else if (comp < 0)
            return binarySearchRec(L, divide(L + M, 2), M, l, m, p, text, suffixArray, lcpArray, seedShape, multiplicity, useLCPTree);
        else
            return binarySearchRec(M, divide(R + M, 2), R, m, r, p, text, suffixArray, lcpArray, seedShape, multiplicity, useLCPTree);

    }

    private final int[] compareResult = new int[2];

    private int[] compare(int m, String p, long tAn, IndexText text, CyclicSeedShape seedShape) {

        int r1, r2;
        do {

            // checking if end of suffix is reached
            int i2 = text.readPosition(tAn + m);
            if (i2 == 127) {
//                int[] result = {m, -1};
//                return result;
                compareResult[0] = m;
                compareResult[1] = -1;
                return compareResult;
            }

            // preparing comparison of m-th pattern character with m-th text character
            if (seedShape.usePosition(m)) {
                r1 = Alphabet.getIndex(Alphabet.reduceCharacter(p.charAt(m)));
                r2 = Alphabet.reducePosition(i2);
            } else { // don't-care-position, continue comparison in any case
                r1 = r2 = 0;
            }
            if (r1 == r2)
                m++;

        } while (r1 == r2 && m < p.length());

//        int[] result = {m, Integer.compare(r1, r2)};
//        return result;

        compareResult[0] = m;
        compareResult[1] = Integer.compare(r1, r2);
        return compareResult;

    }

    private long divide(long num, long denum) {
        return num / denum;
    }

}
