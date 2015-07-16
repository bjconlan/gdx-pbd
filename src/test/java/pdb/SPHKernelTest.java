package pdb;

import com.badlogic.gdx.math.Matrix2;
import com.badlogic.gdx.math.Vector3;
import pbd.SPHKernel;
import org.junit.Test;
import static org.junit.Assert.*;

public class SPHKernelTest {
	SPHKernel r1 = new SPHKernel(1f);
	SPHKernel r2 = new SPHKernel(2f);
	Vector3 a = new Vector3(0f, 1f, 2f);
	Vector3 b = new Vector3(0.3f, 0.3f, 0.3f);

	@Test
	public void wTest() {
		assertEquals(r1.w(a), -9.618275f, 0);
		assertEquals(r1.w(b), 0.564596f, 0);
		assertEquals(r2.w(a), -0.0010468911f, 0);
		assertEquals(r2.w(b), 0.22288762f, 0);
	}

	@Test
	public void gradWTest() {
		assertEquals(r1.gradW(a), new Vector3(new float[] {-0.0f, -10.439774f, -20.879547f}));
		assertEquals(r1.gradW(b), new Vector3(new float[] {-2.0356786f, -2.0356786f, -2.0356786f}));
		assertEquals(r2.gradW(a), new Vector3(new float[] {-0.0f, -0.0059497766f, -0.011899553f}));
		assertEquals(r2.gradW(b), new Vector3(new float[] {-0.1748348f, -0.1748348f, -0.1748348f}));
	}

	@Test
	public void wZeroTest() {
		assertEquals(r1.wZero(), 2.546479f, 0);
		assertEquals(r2.wZero(), 0.31830987f, 0);
	}
}
