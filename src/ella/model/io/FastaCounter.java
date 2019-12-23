package ella.model.io;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipException;

import ella.model.aligner.utils.Alphabet;
import ella.model.aligner.utils.Minimizer;
import ella.model.aligner.utils.streams.MyLongStream;

public class FastaCounter {

	public static Object[] run(File file, int sepDepth, Long size, int p, int q) {
		try {

			long PREFIX_NUMBER = cmpNumOfPrefixes(sepDepth);
			long MEMORY_LIMIT = MyParameters.MAX_MEMORY - 1 * (int) Math.pow(10, 9) - PREFIX_NUMBER * 4;
			long LETTER_LIMIT = (2 * (long) Integer.MAX_VALUE) - 100000;

			if (size != null)
				MEMORY_LIMIT = Math.min(MEMORY_LIMIT, size);

			InputStream is;
			try {
				is = new BufferedInputStream(new GZIPInputStream(new FileInputStream(file)));
			} catch (ZipException e) {
				is = new BufferedInputStream(new FileInputStream(file));
			}

			try {

				// initializing fields for minimizer counting
				Minimizer minimizer = new Minimizer(p, q);
				String pMer = "";
				int seqPos = 0, lastMin = -1;

				// initializing relevant file info
				long totalSeqCounter = 0, totalAACounter = 0;
				MyLongStream batchBounderies = new MyLongStream();
				batchBounderies.add(0);

				// scanning entire file
				int readChars = 0;
				byte[] buffer = new byte[1024];
				long seqCounter = 0, aaCounter = 0, sizeCounter = 200, headerCounter = 0, filePointer = 0;
				boolean empty = true, readingSequence = true, indexHeader = false;
				while ((readChars = is.read(buffer)) != -1) {
					empty = false;
					for (int i = 0; i < readChars; ++i) {
						char c = (char) buffer[i];
						filePointer++;
						if (c == '>') {
							headerCounter += 2;
							sizeCounter += 8;
							seqCounter++;
							readingSequence = false;
							indexHeader = true;
						} else if (readingSequence) {
							if (c != '\n') {
								if (p < q) {
									pMer = extendPMer(pMer, c, p);
									if (pMer.length() == p) {
										minimizer.addPMer(seqPos, pMer);
										int min = (int) minimizer.getMinimium();
										if (seqPos >= q - 1 && min != lastMin) {
											lastMin = min;
											sizeCounter += 7;
										}
									}
								} else
									sizeCounter += 7;
								aaCounter++;
								seqPos++;
							} else {
								pMer = "";
								seqPos = 0;
								lastMin = -1;
								minimizer = new Minimizer(p, q);
							}
						} else if (c == '\n') {
							readingSequence = true;
							indexHeader = false;
						} else if (c == ' ')
							indexHeader = false;
						else if (indexHeader) {
							headerCounter++;
							sizeCounter++;
						}

						if (aaCounter + headerCounter > LETTER_LIMIT || sizeCounter > MEMORY_LIMIT) {
							totalSeqCounter += seqCounter;
							totalAACounter += aaCounter;
							seqCounter = 0;
							aaCounter = 0;
							sizeCounter = 0;
							headerCounter = 0;
							batchBounderies.add(filePointer);
						}

					}
				}

				// finalizing counter and segments
				totalSeqCounter += seqCounter;
				totalAACounter += aaCounter;
				batchBounderies.add(filePointer);
				totalSeqCounter = (totalSeqCounter == 0 && !empty) ? 1 : totalSeqCounter;

				// returning file information
				Object[] result = { totalSeqCounter, totalAACounter, batchBounderies.toArray() };
				return result;

			} finally {
				is.close();
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	private static String extendPMer(String pMer, char c, int p) {
		pMer += c;
		return pMer.length() > p ? pMer.substring(1) : pMer;
	}

	private static int cmpNumOfPrefixes(int sepDepth) {
		int n = Alphabet.getReducedAminoacids().length();
		for (int i = 0; i < sepDepth - 1; i++)
			n *= Alphabet.getReducedAminoacids().length();
		return n;
	}

}
