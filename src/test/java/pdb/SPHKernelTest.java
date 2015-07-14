package pdb;

import com.badlogic.gdx.math.Vector3;
import org.junit.Test;
import pbd.SPHKernel;

public class SPHKernelTest {
	SPHKernel r1 = new SPHKernel(1f);
	SPHKernel r2 = new SPHKernel(2f);
	Vector3 a = new Vector3(0f, 1f, 2f);
	Vector3 b = new Vector3(0.3f, 0.3f, 0.3f);

	@Test
	public void wTest() {
		float r1aResult = r1.w(a);
		float r1bResult = r1.w(b);
		float r2aResult = r2.w(a);
		float r2bResult = r2.w(b);
	}

	@Test
	public void gradWTest() {
		Vector3 r1aResult = r1.gradW(a);
		Vector3 r1bResult = r1.gradW(b);
		Vector3 r2aResult = r2.gradW(a);
		Vector3 r2bResult = r2.gradW(b);

		System.out.println(r1.gradW(a));
	}

	@Test
	public void wZeroTest() {
		float r1Result = r1.wZero();
		float r2Result = r2.wZero();
	}
}
