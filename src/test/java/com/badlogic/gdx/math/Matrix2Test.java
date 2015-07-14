package com.badlogic.gdx.math;


import org.junit.Test;
import static org.junit.Assert.*;

public class Matrix2Test {
	Matrix2 a = new Matrix2(new float[] {1f, 2f, 3f, 4f});

	@Test
	public void idtTest() {
		Matrix2 identity =  new Matrix2(new float[] {1, 0, 0, 1});
		assertEquals(new Matrix2(), identity);
		assertEquals(Matrix2.ZERO.idt(), identity);
	}

	@Test
	public void setMatrix2() {
		Matrix2 mtxZero = Matrix2.ZERO;
		assertEquals(mtxZero.set(new Matrix2(a)), a);
	}

	@Test
	public void detTest() {
		assertEquals(a.det(), -2f, 0);
	}

	@Test
	public void invTest() {
		assertEquals(new Matrix2(a).inv(), new Matrix2(new float[] {-2f, 1.5f, 2f, -0.75f}));
	}

	@Test
	public void transposeTest() {
		assertEquals(new Matrix2(a).transpose(), new Matrix2(new float[] {1f, 3f, 2f, 4f}));
	}
}
