package ella.model.aligner.aligning;

import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicLong;

import ella.model.aligner.indexing.IndexText;
import ella.model.aligner.indexing.SuffixArrayTable;
import ella.model.aligner.utils.CyclicSeedShape;
import ella.model.aligner.utils.bigArrays.UnsignedByteArray;
import ella.model.aligner.utils.bigArrays.UnsignedIntArray;
import ella.model.aligner.utils.suffixArray.EnhancedSuffixArray;
import ella.model.aligner.utils.suffixArray.LCPTreeArray;
import ella.model.aligner.utils.wrapper.ByteArrayWrapper;
import ella.model.io.MyParameters;

public class IndexReader {

	private AtomicLong runtime = new AtomicLong(0);

	private final static long id = 23062016;
	private final static int build = 110;
	private final static int version = 1;

	private IndexText indexText;
	private long numOfSequences, numOfLetters, avgSequenceLength;
	private SuffixArrayTable saTable;
	private ArrayList<SuffixArrayReader> saReaders;

	private File edbFile;
	private FileChannel fC;

	public IndexReader(File edbFile) {
		try {
			this.edbFile = edbFile;
			this.fC = FileChannel.open(edbFile.toPath());
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public boolean run() {

		try {

			if (fC.position() >= Files.size(edbFile.toPath()))
				return false;

			readInHeader();
			readInText();
			setupSuffixArrayReaders();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return true;

	}

	private long readInHeader() throws IOException {

		ByteBuffer buffer = ByteBuffer.allocate(16);
		read(fC, buffer);
		buffer.rewind();
		long magicNumber = buffer.getLong();
		if (magicNumber != id)
			throw new IllegalArgumentException("Not an appropriate index file: " + edbFile.getAbsolutePath());
		int build = buffer.getInt();
		int vers = buffer.getInt();
		if (vers != version)
			System.err.println("WARNING - program version " + version + " does not match index version " + vers + "!");

		return 16;
	}

	private long readInText() throws IOException {

		// parsing text parameters
		ByteBuffer buffer = ByteBuffer.allocate(32);
		read(fC, buffer);
		buffer.rewind();
		MyParameters.TOTAL_BATCHES = buffer.getInt();
		numOfSequences = buffer.getLong();
		numOfLetters = buffer.getLong();
		avgSequenceLength = buffer.getInt();
		int textLength1 = buffer.getInt();
		int textLength2 = buffer.getInt();

		buffer = ByteBuffer.wrap(new byte[textLength1]);
		read(fC, buffer); // fC.read(buffer);
		byte[] text1 = buffer.array();
		buffer = ByteBuffer.wrap(new byte[textLength2]);
		read(fC, buffer);
		byte[] text2 = buffer.array();
		buffer = ByteBuffer.wrap(new byte[8 * (int) numOfSequences]);
		read(fC, buffer);
		byte[] locations = buffer.array();

		// setting up index text
		UnsignedByteArray text = new UnsignedByteArray(text1, text2);
		indexText = new IndexText(text, locations, numOfSequences, numOfLetters);

		return 32L + (long) text1.length + (long) text2.length + (long) locations.length;

	}

	private void read(FileChannel channel, ByteBuffer buffer) {
		try {
			channel.read(buffer);
			while (buffer.hasRemaining())
				channel.read(buffer);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public IndexText getIndexText() {
		return indexText;
	}

	private void setupSuffixArrayReaders() throws IOException {

		// parsing header
		ByteBuffer buffer = ByteBuffer.allocate(3 * 4 + 8);
		read(fC, buffer);
		buffer.rewind();
		int p = buffer.getInt();
		int q = buffer.getInt();
		int sepDepth = buffer.getInt();
		long offset = buffer.getLong();

		saTable = new SuffixArrayTable(indexText, sepDepth, p, q);

		// computing boundaries for suffix array threads
		fC.position(offset);

		saReaders = new ArrayList<SuffixArrayReader>();
		for (int k = 0; k < 2; k++) {

			buffer = ByteBuffer.allocate(4);
			read(fC, buffer);
			buffer.rewind();
			int n = buffer.getInt();

			// reading-in suffix array pointers
			buffer = ByteBuffer.allocate(n * 8);
			read(fC, buffer);
			buffer.rewind();
			ArrayList<Long> saPointers = new ArrayList<Long>(n);
			for (int i = 0; i < n; i++) {
				long pointer = buffer.getLong();
				saPointers.add(pointer);
			}

			// creating suffix array readers
			for (int i = 0; i < saPointers.size(); i++) {
				long l = saPointers.get(i);
				long r = i < saPointers.size() - 1 ? saPointers.get(i + 1) : offset;
				long[] boundary = { l, r };
				saReaders.add(new SuffixArrayReader(sepDepth, boundary, edbFile));
			}

		}

	}

	public ArrayList<SuffixArrayReader> getSaReaders() {
		return saReaders;
	}

	public void freeMemory() {
		for (SuffixArrayReader saReader : saReaders)
			saReader.freeMemory();
		saReaders = null;
	}

	public class SuffixArrayReader {

		private FileChannel saChannel;
		private int sepDepth;
		private long[] boundary;
		private File saFile;

		public SuffixArrayReader(int sepDepth, long[] boundary, File saFile) {
			this.sepDepth = sepDepth;
			this.boundary = boundary;
			this.saFile = saFile;
		}

		public void freeMemory() {
			boundary = null;
		}

		public void run() {

			long time = System.currentTimeMillis();

			try {

				saChannel = FileChannel.open(saFile.toPath());
				long readBytes = boundary[0];
				saChannel.position(boundary[0]);

				try {

					// reading-in prefix
					byte[] bytes = new byte[sepDepth];
					ByteBuffer buffer = ByteBuffer.wrap(bytes);
					read(saChannel, buffer);
					ByteArrayWrapper prefix = new ByteArrayWrapper(buffer.array());
					readBytes += sepDepth;

					// reading in all enhanced suffix arrays corresponding to this prefix
					// buffer.get();
					readByte();
					readBytes += 1;
					while (readBytes < boundary[1]) {

						// parsing seed shape
						StringBuilder builder = new StringBuilder();
						int v;
						while ((v = (int) readByte()) != -96)
							builder.append((char) v);
						CyclicSeedShape seedShape = new CyclicSeedShape(builder.toString());
						readBytes += seedShape.getLength() + 1;

						// parsing suffix array
						byte[] byteSuffixArray1 = new byte[readInt() * 4];
						buffer = ByteBuffer.wrap(byteSuffixArray1);
						read(saChannel, buffer);
						int[] sa1 = byteToIntArray(buffer.array());
						byte[] byteSuffixArray2 = new byte[readInt() * 4];
						buffer = ByteBuffer.wrap(byteSuffixArray2);
						read(saChannel, buffer);
						int[] sa2 = byteToIntArray(buffer.array());
						readInt();
						UnsignedIntArray suffixArray = new UnsignedIntArray(sa1, sa2);
						readBytes += 4 + 4 + suffixArray.size() * 4 + 4;

						// parsing lcp array
						byte[] lcp1 = new byte[readInt()];
						buffer = ByteBuffer.wrap(lcp1);
						read(saChannel, buffer);
						byte[] lcp2 = new byte[readInt()];
						buffer = ByteBuffer.wrap(lcp2);
						read(saChannel, buffer);
						readByte();
						UnsignedByteArray lcpArray = new UnsignedByteArray(lcp1, lcp2);
						LCPTreeArray lcp = new LCPTreeArray(lcpArray, prefix.getAAPrefix());
						readBytes += 4 + 4 + lcpArray.size() + 1;

						// parsing bucket table
						// long n = readLong();
						// int depth = readInt();
						// byte[] byteBuck1 = new byte[readInt() * 4];
						// buffer = ByteBuffer.wrap(byteBuck1);
						// saChannel.read(buffer);
						// int[] buck1 = byteToIntArray(buffer.array());
						// byte[] byteBuck2 = new byte[readInt() * 4];
						// buffer = ByteBuffer.wrap(byteBuck2);
						// saChannel.read(buffer);
						// int[] buck2 = byteToIntArray(buffer.array());
						// readByte();
						// UnsignedIntArray buckets = new UnsignedIntArray(buck1, buck2);
						// BucketTable bucketTable = new BucketTable(n, depth, buckets, seedShape);
						// readBytes += 8 + 4 + 4 + 4 + buckets.size() * 4 + 1;

						// setting up enhanced suffix array
						EnhancedSuffixArray esa = new EnhancedSuffixArray(prefix, seedShape, suffixArray, lcp, null);
						saTable.addSuffixArray(esa);

					}

					buffer = null;

				} finally {
					saChannel.close();
				}
			} catch (Exception e) {
				e.printStackTrace();
			}

			runtime.getAndAdd(System.currentTimeMillis() - time);

		}

		private int readByte() throws IOException {
			ByteBuffer buffer = ByteBuffer.allocate(1);
			read(saChannel, buffer);
			buffer.rewind();
			return buffer.get();
		}

		private int readInt() throws IOException {
			ByteBuffer buffer = ByteBuffer.allocate(4);
			read(saChannel, buffer);
			buffer.rewind();
			return buffer.getInt();
		}

		private long readLong() throws IOException {
			ByteBuffer buffer = ByteBuffer.allocate(8);
			read(saChannel, buffer);
			buffer.rewind();
			return buffer.getLong();
		}

	}

	public int[] byteToIntArray(byte buf[]) {
		int intArr[] = new int[buf.length / 4];
		int offset = 0;
		for (int i = 0; i < intArr.length; i++) {
			intArr[i] = (buf[3 + offset] & 0xFF) | ((buf[2 + offset] & 0xFF) << 8) | ((buf[1 + offset] & 0xFF) << 16)
					| ((buf[0 + offset] & 0xFF) << 24);
			offset += 4;
		}
		return intArr;
	}

	public void close() {
		try {
			fC.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public SuffixArrayTable getSaTable() {
		return saTable;
	}

	public long getNumOfSequences() {
		return numOfSequences;
	}

	public long getNumOfLetters() {
		return numOfLetters;
	}

	public long getAvgSequenceLength() {
		return avgSequenceLength;
	}

	public long getRuntime() {
		return runtime.get();
	}

}
