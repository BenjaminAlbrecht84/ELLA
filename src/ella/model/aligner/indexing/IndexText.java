package ella.model.aligner.indexing;

import java.nio.ByteBuffer;

import ella.model.aligner.utils.Alphabet;
import ella.model.aligner.utils.bigArrays.UnsignedByteArray;
import ella.model.aligner.utils.bigArrays.UnsignedIntArray;

public class IndexText {

	private final long numOfSeqs, numOfLetters;
	private int avgSequenceLength;
	private final UnsignedByteArray text;
	private final int[][] locations;

	public IndexText(UnsignedByteArray text, UnsignedIntArray locationInfo, long numOfSeqs, long numOfLetters, int avgSequenceLength) {
		this.text = text;
		this.numOfSeqs = numOfSeqs;
		this.numOfLetters = numOfLetters;
		this.avgSequenceLength = avgSequenceLength;
		locations = parseLocations(locationInfo);
	}

	private int[][] parseLocations(UnsignedIntArray locationInfo) {
		int[][] locations = new int[2][(int) (locationInfo.size() / 2)];
		int pos = 0;
		for (long i = 0; i < locationInfo.size(); i++) {
			int value = (int) locationInfo.get(i);
			int j = (i + 1) % 2 == 0 ? 1 : 0;
			locations[j][pos] = value;
			if (j == 1)
				pos++;
		}
		return locations;
	}

	public IndexText(UnsignedByteArray text, byte[] locationInfo, long numOfSeqs, long numOfLetters) {
		this.text = text;
		this.numOfSeqs = numOfSeqs;
		this.numOfLetters = numOfLetters;
		locations = parseLocations(locationInfo);
	}

	private int[][] parseLocations(byte[] locationInfo) {
		int[][] locations = new int[2][locationInfo.length / 8];
		byte[] b = new byte[4];
		int pos = 0;
		for (int i = 0; i < locationInfo.length; i++) {
			b[i % 4] = locationInfo[i];
			if ((i + 1) % 4 == 0) {
				int value = ByteBuffer.wrap(b).getInt();
				int j = (i + 1) % 8 == 0 ? 1 : 0;
				locations[j][pos] = value;
				if (j == 1)
					pos++;
			}
		}
		return locations;
	}

	public int getProteinLocationIndex(long pos) {
		int i = binarySearch(locations[0], pos);
		i = i < 0 ? -i - 2 : i;
		return i;
	}

	public Object[] getProteinInfo(long pos) {
		int i = binarySearch(locations[0], pos);
		i = i < 0 ? -i - 2 : i;
		Object[] info = { Integer.toUnsignedLong(locations[0][i]), locations[1][i] };
		return info;
	}

	public long getProteinStart(long pos) {
		int i = binarySearch(locations[0], pos);
		i = i < 0 ? -i - 2 : i;
		return Integer.toUnsignedLong(locations[0][i]);
	}

	public int getProteinLength(long pos) {
		int i = binarySearch(locations[0], pos);
		i = i < 0 ? -i - 2 : i;
		return locations[1][i];
	}

	private static int binarySearch(int[] a, long key) {
		return binarySearchRec(a, 0, a.length, key);
	}

	private static int binarySearchRec(int[] a, int fromIndex, int toIndex, long key) {
		int low = fromIndex;
		int high = toIndex - 1;
		while (low <= high) {
			int mid = (low + high) >>> 1;
			long midVal = Integer.toUnsignedLong(a[mid]);
			if (midVal < key)
				low = mid + 1;
			else if (midVal > key)
				high = mid - 1;
			else
				return mid; // key found
		}
		return -(low + 1); // key not found.
	}

	public int getProteinLengthFromIndex(int i) {
		return locations[1][i];
	}

	public String getProteinSequencen(long pos) {
		long i = getProteinStart(pos);
		StringBuffer buf = new StringBuffer();
		while (readPosition(i) != 127)
			buf.append(readAA(i++));
		return buf.toString();
	}

	public String getProteinAccession(long pos) {
		long start = getProteinStart(pos);
		int length = getProteinLength(pos);
		long i = start + length + 1;
		StringBuffer buf = new StringBuffer();
		while (readPosition(i) != 127)
			buf.append(readCharacter(i++));
		return buf.toString();
	}

	public String getProteinAccessionFromIndex(int i) {
		StringBuffer buf = new StringBuffer();
		long pos = Integer.toUnsignedLong(locations[0][i]) + Integer.toUnsignedLong(locations[1][i]) + 1;
		while (readPosition(pos) != 127)
			buf.append(readCharacter(pos++));
		return buf.toString();
	}

	public int readPosition(long pos) {
		return text.get(pos);
	}

	public char readCharacter(long pos) {
		char c = (char) readPosition(pos);
		return c;
	}

	public char readAA(long pos) {
		char aa = Alphabet.getCharacter(readPosition(pos));
		return aa;
	}

	public UnsignedByteArray getText() {
		return text;
	}

	public long length() {
		return text.size();
	}

	public long getByteSize() {
		return text.size() + locations[0].length * 4 * 2;
	}

	public long getNumOfSeqs() {
		return numOfSeqs;
	}

	public long getNumOfLetters() {
		return numOfLetters;
	}

	public int getAvgSequenceLength() {
		return avgSequenceLength;
	}

	public int[][] getLocations() {
		return locations;
	}

	public int[] getLocationsAsArray() {
		int[] a = new int[locations[0].length * 2];
		int pos = 0;
		for (int j = 0; j < locations[0].length; j++) {
			for (int i = 0; i < locations.length; i++) {
				a[pos++] = locations[i][j];
			}
		}
		return a;
	}

	public String toFASTA() {
		StringBuffer acc = new StringBuffer();
		StringBuffer seq = new StringBuffer();
		StringBuffer buf = new StringBuffer();
		for (int i = 0; i < text.size(); i++) {
			int v = readPosition(i);
			if (v < Alphabet.getAaString().length())
				seq.append(readAA(i));
			else if (v != 127)
				acc.append(readCharacter(i));
			else if (acc.length() > 0) {
				buf.append(">" + acc.toString() + "\n");
				buf.append(seq.toString() + "\n");
				acc = new StringBuffer();
				seq = new StringBuffer();
			}
		}
		return buf.toString();
	}

	public String toString() {
		StringBuffer buf = new StringBuffer();
		for (int i = 0; i < text.size(); i++) {
			int v = readPosition(i);
			if (v < Alphabet.getAaString().length()) {
				buf.append(readAA(i));
			} else if (v == 127) {
				buf.append("$");
			} else
				buf.append(this.readCharacter(i));
		}
		for (int j = 0; j < locations[0].length; j++) {
			for (int i = 0; i < locations.length; i++) {
				buf.append("[" + locations[i][j] + "]\n");
			}
		}
		return buf.toString();
	}

}
