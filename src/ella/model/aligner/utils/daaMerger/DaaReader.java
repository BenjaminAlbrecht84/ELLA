/*
 * Copyright 2017 Benjamin Albrecht
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

package ella.model.aligner.utils.daaMerger;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicLong;

import ella.model.aligner.utils.SparseString;

public class DaaReader {

	private File daaFile;
	private boolean verbose = false;

	private DaaHeader header;
	private CountDownLatch latch;
	private ConcurrentHashMap<String, ReadHits> readMap;

	private AtomicLong allParsedRecords;;
	private int last_p;
	private long numOfRecords;

	public DaaReader(File daaFile, boolean verbose) {
		this.verbose = verbose;
		this.daaFile = daaFile;
		header = new DaaHeader(daaFile);
		header.loadAllReferences();
		if (verbose)
			header.print();
	}

	public ConcurrentHashMap<String, ReadHits> parseAllHits(int cores) {

		System.out.println("STEP_3>Parsing DIAMOND output...");
		long time = System.currentTimeMillis();

		readMap = new ConcurrentHashMap<>();

		System.out.println("OUTPUT>Parsing " + header.getNumberOfQueryRecords() + " query records...");

		last_p = 0;
		allParsedRecords = new AtomicLong(0);
		numOfRecords = header.getNumberOfQueryRecords();

		int chunk = (int) Math.ceil((double) header.getNumberOfQueryRecords() / (double) cores);
		Vector<Thread> allParser = new Vector<Thread>();
		for (int l = 0; l < header.getNumberOfQueryRecords(); l += chunk) {
			int r = (l + chunk) < header.getNumberOfQueryRecords() ? l + chunk : (int) header.getNumberOfQueryRecords();
			int[] bounds = { l, r };
			DaaParser parser = new DaaParser(bounds);
			allParser.add(parser);
		}

		latch = new CountDownLatch(allParser.size());
		ExecutorService executor = Executors.newFixedThreadPool(cores);
		for (Thread thread : allParser)
			executor.execute(thread);

		// awaiting termination
		try {
			latch.await();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		executor.shutdown();

		long runtime = (System.currentTimeMillis() - time) / 1000;
		System.out.println("OUTPUT>" + 100 + "% (" + numOfRecords + "/" + numOfRecords + ") of all records parsed.[" + runtime + "s]\n");

		return readMap;

	}

	private synchronized void reportProgress(int p) {
		p = ((int) Math.floor((double) p / 10.)) * 10;
		if (p != 100 && p != last_p && p % 1 == 0) {
			System.out.println("OUTPUT>" + p + "% (" + allParsedRecords + "/" + numOfRecords + ") of all records parsed.");
			last_p = p;
		}
	}

	public QueryIterator queryIterator() throws IOException {
		return new QueryIterator();
	}

	public class QueryIterator implements Iterator<QueryHits>, Closeable {

		final private RandomAccessFile raf;
		private long queryPointer = 0;
		private final long totalQueries;

		public QueryIterator() throws IOException {
			totalQueries = header.getNumberOfQueryRecords();
			raf = new RandomAccessFile(daaFile, "r");
			raf.seek(header.getLocationOfBlockInFile(header.getAlignmentsBlockIndex()));
		}

		public long size() {
			return totalQueries;
		}

		public boolean hasNext() {
			return queryPointer < totalQueries;
		}

		public QueryHits next() {
			try {
				if (queryPointer >= totalQueries)
					return null;
				queryPointer++;

				QueryHits queryHits = new QueryHits();
				DaaHit hit = new DaaHit();

				long filePointer = raf.getFilePointer();
				ByteBuffer buffer = ByteBuffer.allocate(4);
				buffer.order(ByteOrder.LITTLE_ENDIAN);
				raf.read(buffer.array());
				int alloc = buffer.getInt();

				ByteBuffer hitBuffer = ByteBuffer.allocate(alloc);
				hitBuffer.order(ByteOrder.LITTLE_ENDIAN);
				raf.read(hitBuffer.array());

				// parsing query properties
				hit.parseQueryProperties(filePointer, hitBuffer, false, true);
				queryHits.setQueryLength(hit.getTotalQueryLength());
				queryHits.setQueryName(hit.getQueryName());

				while (hitBuffer.position() < hitBuffer.capacity()) {

					int accessPoint = hitBuffer.position();

					// parsing match properties
					hit.parseHitProperties(header, hitBuffer, false);

					int refStart = hit.getRefStart() + 1;
					int bitScore = hit.getBitScore();
					int rawScore = hit.getRawScore();
					int queryStart = hit.getQueryStart() + 1;
					int refLength = hit.getTotalRefLength();
					int queryLength = hit.getQueryLength();
					int refCover = hit.getRefLength();
					int totalQueryLength = hit.getTotalQueryLength();
					int frame = hit.getFrame();
					SparseString subject = new SparseString(hit.getReferenceName().split(" ")[0]);
					int subjectID = hit.getSubjectID();
					String query = hit.getQueryName();
					String refName = new String(getDaaHeader().getReferenceName(subjectID));

					// initializing hit
					AlignmentHit h = new AlignmentHit(refStart, bitScore, rawScore, queryStart, refLength, totalQueryLength, subjectID,
							hit.getEditByteOperations(), frame, query, refName, filePointer, accessPoint);
					h.setFrame(frame);

					// storing hit
					queryHits.add(h, subject, frame);
				}
				return queryHits;
			} catch (IOException ex) {
				throw new RuntimeException(ex);
			}
		}

		@Override
		public void close() throws IOException {
			raf.close();
		}
	}

	private synchronized void addHits(HashMap<String, ReadHits> localReadMap) {
		for (String read_id : localReadMap.keySet()) {
			if (!readMap.containsKey(read_id))
				readMap.put(read_id, new ReadHits());
			ReadHits hits = localReadMap.get(read_id);
			for (SparseString subject : hits.getHitMap().keySet()) {
				for (int frame : hits.getHitMap().get(subject).keySet()) {
					for (Hit h : hits.getHitMap().get(subject).get(frame))
						readMap.get(read_id).add(h, subject, frame);
				}
			}
		}
	}

	public Object[] parseDAAHitByIndex(int index, Integer lastIndex, Long lastFilePointer) {

		ArrayList<DaaHit> daaHits = new ArrayList<DaaHit>();

		try {

			RandomAccessFile raf = new RandomAccessFile(daaFile, "r");
			raf.seek(header.getLocationOfBlockInFile(header.getAlignmentsBlockIndex()));

			int start = lastIndex != null ? lastIndex : 0;
			if (lastFilePointer != null)
				raf.seek(lastFilePointer);
			int lastIteration = 0;
			long filePointer = raf.getFilePointer();
			for (int i = start; i < header.getNumberOfQueryRecords(); i++) {

				lastIteration = i;
				filePointer = raf.getFilePointer();

				ByteBuffer buffer = ByteBuffer.allocate(4);
				buffer.order(ByteOrder.LITTLE_ENDIAN);
				raf.read(buffer.array());
				int alloc = buffer.getInt();

				ByteBuffer hitBuffer = ByteBuffer.allocate(alloc);
				hitBuffer.order(ByteOrder.LITTLE_ENDIAN);
				raf.read(hitBuffer.array());

				if (i == index) {
					DaaHit hit = new DaaHit();
					hit.parseQueryProperties(filePointer, hitBuffer, false, true);
					while (hitBuffer.position() < hitBuffer.capacity()) {					
						DaaHit h = new DaaHit();
						h.copyQueryProperties(hit);
						h.parseHitProperties(header, hitBuffer, false);
						daaHits.add(h);
					}
					break;
				}

			}

			raf.close();
			Object[] result = { daaHits, lastIteration, filePointer };
			return result;

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		Object[] result = { daaHits, null, null };
		return result;

	}

	public class DaaParser extends Thread {

		private int[] bounds;

		public DaaParser(int[] queryRecordsBounds) {
			this.bounds = queryRecordsBounds;
		}

		public void run() {

			HashMap<String, ReadHits> localReadMap = new HashMap<String, ReadHits>();

			try {

				RandomAccessFile raf = new RandomAccessFile(daaFile, "r");

				try {

					raf.seek(header.getLocationOfBlockInFile(header.getAlignmentsBlockIndex()));
					for (int i = 0; i < header.getNumberOfQueryRecords(); i++) {

						DaaHit hit = new DaaHit();

						long filePointer = raf.getFilePointer();
						ByteBuffer buffer = ByteBuffer.allocate(4);
						buffer.order(ByteOrder.LITTLE_ENDIAN);
						raf.read(buffer.array());
						int alloc = buffer.getInt();

						ByteBuffer hitBuffer = ByteBuffer.allocate(alloc);
						hitBuffer.order(ByteOrder.LITTLE_ENDIAN);
						raf.read(hitBuffer.array());

						if (i >= bounds[0] && i < bounds[1]) {

							// parsing query properties
							hit.parseQueryProperties(filePointer, hitBuffer, false, false);

							while (hitBuffer.position() < hitBuffer.capacity()) {

								// parsing match properties
								hit.parseHitProperties(header, hitBuffer, false);

								int refStart = hit.getRefStart() + 1;
								int bitScore = hit.getBitScore();
								int rawScore = hit.getRawScore();
								int queryStart = hit.getQueryStart() + 1;
								int refLength = hit.getTotalRefLength();
								int queryLength = hit.getQueryLength() / 3;
								int refCover = hit.getRefLength();

								int frame = hit.getFrame();
								SparseString subject = new SparseString(hit.getReferenceName().split(" ")[0]);
								int subjectID = hit.getSubjectID();

								// parsing query name
								String query = hit.getQueryName();

								// initializing hit
								Hit h = new Hit(refStart, bitScore, rawScore, queryStart, refLength, queryLength,
										subjectID, refCover, hit.getEditByteOperations());
								h.setFrame(frame);

								// storing hit
								if (!localReadMap.containsKey(query))
									localReadMap.put(query, new ReadHits());
								localReadMap.get(query).add(h, subject, frame);

							}

							if (i != 0 && i % 1000 == 0) {
								int p = (int) Math.round(((double) allParsedRecords.addAndGet(1000) / (double) numOfRecords) * 100.);
								reportProgress(p);
							}

						}

						if (i >= bounds[1])
							break;

					}

				} finally {
					raf.close();
				}

			} catch (Exception e) {
				e.printStackTrace();
			}

			addHits(localReadMap);
			latch.countDown();

		}

		private String[] mySplit(String s, char c) {
			List<String> words = new ArrayList<String>();
			int pos = 0, end;
			while ((end = s.indexOf(c, pos)) >= 0) {
				words.add(s.substring(pos, end));
				pos = end + 1;
			}
			if (pos < s.length())
				words.add(s.substring(pos, s.length()));
			String[] entries = words.toArray(new String[words.size()]);
			return entries;

		}

	}

	public DaaHit parseHit(RandomAccessFile raf, long filePointer, int accessPoint) {

		DaaHit hit = null;

		try {

			raf.seek(filePointer);

			ByteBuffer buffer = ByteBuffer.allocate(4);
			buffer.order(ByteOrder.LITTLE_ENDIAN);
			raf.read(buffer.array());
			int alloc = buffer.getInt();

			ByteBuffer hitBuffer = ByteBuffer.allocate(alloc);
			hitBuffer.order(ByteOrder.LITTLE_ENDIAN);
			raf.read(hitBuffer.array());

			// parsing query properties
			hit = new DaaHit();
			hit.parseQueryProperties(filePointer, hitBuffer, true, false);

			hitBuffer.position(accessPoint);

			// parsing match properties
			hit.parseHitProperties(header, hitBuffer, true);

		} catch (Exception e) {
			e.printStackTrace();
		}

		return hit;

	}

	public File getDaaFile() {
		return daaFile;
	}

	public DaaHeader getDaaHeader() {
		return header;
	}

}
