package ella.model.aligner.utils;

import java.util.ArrayList;
import java.util.HashMap;

public class AlignmentCompressor {

	private static HashMap<Character, Integer> aaToIndex = new HashMap<Character, Integer>();
	static {
		String aaString = Alphabet.getAaString();
		for (int i = 0; i < aaString.length(); i++)
			aaToIndex.put(aaString.charAt(i), i);
	}

	public static ArrayList<Byte> run(String[] ali) {

		ArrayList<Byte> editOps = new ArrayList<Byte>();
		char lastType = '-';
		int num = 0;
		for (int i = 0; i < ali[0].length(); i++) {

			char c1 = ali[0].charAt(i);
			char c2 = ali[1].charAt(i);
			char type = getEditType(c1, c2);
			if (type == 'M' || type == 'I') {
				if (type != lastType && num != 0) {
					editOps.addAll(getEditOperation(lastType, num));
					num = 0;
				}
				num++;
			} else {
				if (num != 0 && lastType != '-') {
					editOps.addAll(getEditOperation(lastType, num));
					num = 0;
				}
				editOps.addAll(getEditOperation(type, c2));
			}
			lastType = type;

		}
		if (lastType == 'M' || lastType == 'I')
			editOps.addAll(getEditOperation(lastType, num));

		return editOps;
	}

	private static ArrayList<Byte> getEditOperation(char type, int total) {
		ArrayList<Byte> opVec = new ArrayList<Byte>();
		while (total > 0) {
			int num = total > 63 ? 63 : total;
			total -= num;
			byte op = 0;
			if (type == 'I')
				op |= 1 << 6;
			op |= num;
			opVec.add(op);	
		}
		return opVec;
	}
	
	public static ArrayList<Byte> cmpFrameShiftOperations(int queryDiff) {
		ArrayList<Byte> ops = new ArrayList<Byte>();
		char slash = queryDiff > 0 ? '\\' : '/';
		int length = Math.abs(queryDiff);
		for (int i = 0; i < length; i++) {
			byte op = 0;
			op |= 3 << 6;
			op |= Alphabet.getIndex(slash);
			ops.add(op);
		}
		return ops;
	}

	private static ArrayList<Byte> getEditOperation(char type, char c) {
		byte op = 0;
		if (type == 'D')
			op |= 2 << 6;
		else
			op |= 3 << 6;
		op |= aaToIndex.get(c);
		ArrayList<Byte> opVec = new ArrayList<Byte>();
		opVec.add(op);	
		return opVec;
	}

	private static char getEditType(char c1, char c2) {
		if (c1 == '-')
			return 'D';
		if (c2 == '-')
			return 'I';
		if (c1 != c2)
			return 'S';
		return 'M';
	}

}
