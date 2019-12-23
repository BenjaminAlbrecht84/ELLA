package ella.model.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipException;

import ella.model.aligner.utils.DNAPacker;
import ella.model.aligner.utils.SparseString;

public class FastaReader {

	public static ArrayList<Object[]> read(File fastAQFile) {

		ArrayList<Object[]> readInfo = new ArrayList<Object[]>();
		HashSet<SparseString> readNames = new HashSet<SparseString>();
		try {

			BufferedReader buf;
			try {
				buf = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fastAQFile))));
			} catch (ZipException e) {
				buf = new BufferedReader(new FileReader(fastAQFile));
			}

			String line, id = "";
			boolean readSequence = false;
			StringBuilder seq = new StringBuilder("");
			while ((line = buf.readLine()) != null) {
				if (line.startsWith("@") || line.startsWith(">")) {
					if (seq.length() != 0 && !id.isEmpty()) {						
						Object[] o = { new SparseString(id), DNAPacker.packSequence(seq.toString()), seq.length() };
						if (checkFASTFile(o, readNames))
							readInfo.add(o);
					}
					seq = new StringBuilder("");
					id = line.substring(1).split(" ")[0];
					readSequence = true;
				} else if (line.startsWith("+")) {
					readSequence = false;
				} else if (readSequence) {
					seq.append(line);
				}

			}
			if (seq.length() != 0 && !id.isEmpty()) {
				Object[] o = { new SparseString(id), DNAPacker.packSequence(seq.toString()), seq.length() };
				if (checkFASTFile(o, readNames))
					readInfo.add(o);
			}
			buf.close();

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		readNames = null;

		return readInfo;

	}

	private static boolean checkFASTFile(Object[] o, HashSet<SparseString> readNames) {
		if (readNames.contains(o[0])) {
			System.err.println("WARNING: read " + o[0].toString() + " occurs multiple times in FASTA file.");
			return false;
		}
		readNames.add((SparseString) o[0]);
		return true;
	}

}
