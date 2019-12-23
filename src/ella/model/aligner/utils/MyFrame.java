package ella.model.aligner.utils;

public enum MyFrame {

	FWD_1("+1"), FWD_2("+2"), FWD_3("+3"), REV_1("-1"), REV_2("-2"), REV_3("-3");

	private final String id;
	private final MyStrand strand;

	private MyFrame(String s) {
		id = s;
		strand = Integer.parseInt(id) > 0 ? MyStrand.FWD : MyStrand.REV;
	}

	public String getID() {
		return id;
	}

	public int getNumericalID() {
		return Integer.parseInt(id);
	}

	public MyStrand getStrand() {
		return strand;
	}

}