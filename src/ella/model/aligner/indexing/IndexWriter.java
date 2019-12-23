package ella.model.aligner.indexing;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.file.Files;
import java.util.concurrent.ConcurrentHashMap;

import ella.model.aligner.utils.CyclicSeedShape;
import ella.model.aligner.utils.streams.MyByteOutputStream;
import ella.model.aligner.utils.streams.MyByteStream;
import ella.model.aligner.utils.suffixArray.EnhancedSuffixArray;
import ella.model.aligner.utils.suffixArray.LCPTreeArray;
import ella.model.aligner.utils.wrapper.ByteArrayWrapper;
import ella.model.io.ByteWriter;

public class IndexWriter {

	public final static int MAX_STREAM_CAPACITY = 1000000;

	public static void writeIndexFile2File(IndexText text, File f, long id, int build, int version, boolean firstBatch, int totalBatches)
			throws IOException {

		// initializing offsets
		pointersOut = new MyByteStream();
		textOffset = headerOffset = saOffset = 0;
		batchOffset = firstBatch ? 0 : Files.size(f.toPath());

		long writtenBytes = writeHeader2File(text, f, id, build, version, !firstBatch, totalBatches);
		writtenBytes += writeText2File(text, f, id, build, version);
		textOffset = writtenBytes;

	}

	private static long writeHeader2File(IndexText text, File f, long id, int build, int version, boolean append, int totalBatches)
			throws IOException {

		MyByteStream out = new MyByteStream(40);

		// writing header
		ByteWriter.write(out, id);
		ByteWriter.write(out, build);
		ByteWriter.write(out, version);
		ByteWriter.write(out, totalBatches);
		ByteWriter.write(out, text.getNumOfSeqs());
		ByteWriter.write(out, text.getNumOfLetters());
		ByteWriter.write(out, text.getAvgSequenceLength());

		// writing to file
		writeBytes2File((byte[]) out.toArray().getArrays()[0], f, append);

		return out.size();

	}

	private static long writeText2File(IndexText text, File f, long id, int build, int version) throws IOException {

		MyByteOutputStream out = new MyByteOutputStream(f, MAX_STREAM_CAPACITY, true);

		// writing text
		Object[] textArrays = text.getText().getArrays();
		byte[] text1 = (byte[]) textArrays[0], text2 = (byte[]) textArrays[1];
		ByteWriter.write(out, text1.length);
		ByteWriter.write(out, text2.length);
		ByteWriter.write(out, text1);
		ByteWriter.write(out, text2);

		// writing sequence locations
		ByteWriter.write(out, text.getLocationsAsArray());

		// writing into file
		out.writeToFile();

		return out.getWrittenBytes();

	}

	private static MyByteStream pointersOut = new MyByteStream();
	private static long saOffsetPointer = 0;
	private static long batchOffset = 0, textOffset = 0, headerOffset = 0, saOffset = 0;

	public static void writeSuffixTableHeader2File(int p, int q, int sepDepth, File f) throws IOException {

		// writing header
		MyByteStream headerOut = new MyByteStream(16);
		ByteWriter.write(headerOut, p);
		ByteWriter.write(headerOut, q);
		ByteWriter.write(headerOut, sepDepth);
		saOffsetPointer = batchOffset + textOffset + headerOut.size();
		ByteWriter.write(headerOut, (long) 0);
		headerOffset = headerOut.size();
		writeBytes2File((byte[]) headerOut.toArray().getArrays()[0], f, true);

	}

	public static void writeSuffixArrayTable2File(File outFile, SuffixArrayTable suffixArrayTable, IndexText indexText) throws IOException {

		// writing entire suffix array table to byte stream
		MyByteOutputStream saOut = new MyByteOutputStream(outFile, MAX_STREAM_CAPACITY, true);
		for (ByteArrayWrapper p : suffixArrayTable.getSuffixArrayTable().keySet())
			writeSuffixArrays2File(outFile, p, suffixArrayTable.getSuffixArrayTable().get(p), saOut, indexText);

		// writing stream to file and freeing memory
		saOut.writeToFile();
		saOffset += saOut.getWrittenBytes();
		saOut = null;

	}

	public static void writeSuffixArrays2File(File f, ByteArrayWrapper prefix, ConcurrentHashMap<CyclicSeedShape, EnhancedSuffixArray> saTable,
			MyByteOutputStream saOut, IndexText indexText) throws IOException {

		// writing suffix arrays grouped by prefixes and seed-shapes
		long saPointer = batchOffset + textOffset + headerOffset + saOffset + saOut.getWrittenBytes();
		ByteWriter.write(pointersOut, saPointer);
		ByteWriter.write(saOut, prefix.getData());
		ByteWriter.write(saOut, (byte) -96);

		for (CyclicSeedShape seedShape : saTable.keySet()) {

			// writing seedShape
			ByteWriter.write(saOut, seedShape.toString());
			ByteWriter.write(saOut, (byte) -96);

			// writing suffix array
			EnhancedSuffixArray esa = saTable.get(seedShape);
			int[] esa1 = (int[]) esa.getSuffixArray().getArrays()[0], esa2 = (int[]) esa.getSuffixArray().getArrays()[1];
			ByteWriter.write(saOut, esa1.length);
			ByteWriter.write(saOut, esa1);
			ByteWriter.write(saOut, esa2.length);
			ByteWriter.write(saOut, esa2);
			ByteWriter.write(saOut, (int) -96);

			// writing LCP-array
			LCPTreeArray lcp = esa.getLcpTreeArray();
			byte[] lcp1 = (byte[]) lcp.getLcpTreeArray().getArrays()[0], lcp2 = (byte[]) lcp.getLcpTreeArray().getArrays()[1];
			ByteWriter.write(saOut, lcp1.length);
			ByteWriter.write(saOut, lcp1);
			ByteWriter.write(saOut, lcp2.length);
			ByteWriter.write(saOut, lcp2);
			ByteWriter.write(saOut, (byte) -96);

			// writing bucket-table
			// BucketTable bucketTable = esa.getBucketTable();
			// ByteWriter.write(saOut, bucketTable.getN());
			// ByteWriter.write(saOut, bucketTable.getDepth());
			// int[] buck1 = (int[]) bucketTable.toArray().getArrays()[0], buck2 = (int[]) bucketTable.toArray().getArrays()[1];
			// ByteWriter.write(saOut, buck1.length);
			// ByteWriter.write(saOut, buck1);
			// ByteWriter.write(saOut, buck2.length);
			// ByteWriter.write(saOut, buck2);
			// ByteWriter.write(saOut, (byte) -96);

		}

	}

	public static void writePointers2File(File f) throws IOException {
		RandomAccessFile raf = new RandomAccessFile(f, "rw");
		raf.seek(saOffsetPointer);

		long offset = batchOffset + textOffset + headerOffset + saOffset;
		raf.writeLong(offset);
		raf.close();

		byte[] bytes0 = (byte[]) pointersOut.toArray().getArrays()[0];
		writeInteger2File(bytes0.length / 8, f, true);
		writeBytes2File(bytes0, f, true);

		byte[] bytes1 = (byte[]) pointersOut.toArray().getArrays()[1];
		writeInteger2File(bytes1.length / 8, f, true);
		writeBytes2File(bytes1, f, true);

	}

	private static void writeInteger2File(int value, File f, boolean append) throws IOException {
		byte[] bytes = new byte[4];
		bytes[0] = (byte) ((value & 0xFF000000) >> 24);
		bytes[1] = (byte) ((value & 0x00FF0000) >> 16);
		bytes[2] = (byte) ((value & 0x0000FF00) >> 8);
		bytes[3] = (byte) ((value & 0x000000FF) >> 0);
		writeBytes2File(bytes, f, append);
	}

	private static void writeBytes2File(byte[] bytes, File f, boolean append) throws IOException {
		FileOutputStream fos = new FileOutputStream(f.getAbsolutePath(), append);
		fos.write(bytes);
		fos.close();
	}

}
