package com.badlogic.gdx.math;

import com.badlogic.gdx.utils.GdxRuntimeException;

import java.io.Serializable;
import java.util.Arrays;

public class Matrix2 implements Serializable {
	private static final long serialVersionUID = 1L;
	public static final int M00 = 0;
	public static final int M01 = 1;
	public static final int M10 = 2;
	public static final int M11 = 3;
	public float[] val = new float[4];

	public Matrix2() {
		idt();
	}

	public Matrix2(Matrix2 matrix) {
		set(matrix);
	}

	public Matrix2(float[] values) {
		this.set(values);
	}

	public Matrix2 idt() {
		float[] val = this.val;
		val[M01] = val[M10] = 0;
		val[M00] = val[M11] = 1;
		return this;
	}

	public Matrix2 set(Matrix2 mat) {
		System.arraycopy(mat.val, 0, val, 0, val.length);
		return this;
	}

	public Matrix2 set(float[] values) {
		System.arraycopy(values, 0, val, 0, val.length);
		return this;
	}

	// determinant
	public float det() {
		return val[M00] * val[M11] - val[M10] * val[M01];
	}

	public Matrix2 inv() {
		float det = det();
		if (det != 0) det = 1 / det;
		val[M00] = det * val[M11];
		val[M01] = -det * val[M10];
		val[M10] = -det * val[M11];
		val[M11] = det * val[M01];

		return this;
	}

	public Matrix2 transpose() {
		float v01 = val[M10];
		float v10 = val[M01];
		val[M01] = v01;
		val[M10] = v10;
		return this;
	}

	// used for loops but best to use MXX statics as accessors
	public float get(int col, int row) {
		return val[(col | row << 1) & 3]; // & mask
	}

	// FIXME use kotlin to extend Affine2.set(Matrix2)
	public Affine2 toAffine2() {
		Affine2 r = new Affine2();
		r.m00 = val[M00]; r.m01 = val[M01]; r.m02 = 0;
		r.m10 = val[M10]; r.m11 = val[M11]; r.m02 = 0;
		return r;
	}

	@Override
	public String toString() {
		return "Matrix2{" +
				"val=" + Arrays.toString(val) +
				'}';
	}
}
