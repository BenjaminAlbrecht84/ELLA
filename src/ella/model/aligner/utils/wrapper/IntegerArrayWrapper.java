package ella.model.aligner.utils.wrapper;

import java.util.Arrays;

public class IntegerArrayWrapper {

	private final int[] data;

	public IntegerArrayWrapper(int[] data) {
		if (data == null) {
			throw new NullPointerException();
		}
		this.data = data;
	}

	@Override
	public boolean equals(Object o) {
		if (!(o instanceof IntegerArrayWrapper)) {
			return false;
		}
		return Arrays.equals(data, ((IntegerArrayWrapper) o).data);
	}

	@Override
	public int hashCode() {
		return Arrays.hashCode(data);
	}

	public int[] getData() {
		return data;
	}

}
