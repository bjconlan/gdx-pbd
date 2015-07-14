import com.badlogic.gdx.math.Matrix2;
import org.junit.Test;

public class Test1 {

	@Test
	public void sample() {
		System.out.println((1 | 1 << 1) & 3);
		//System.out.println(1 << 1 & 1 << 1);
		Matrix2 m2 = new Matrix2();
		m2.set(0, 1, 3.0f);
		System.out.println(m2);
	}
}
