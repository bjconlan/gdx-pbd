package pdb;

import com.badlogic.gdx.math.Matrix3;
import com.badlogic.gdx.math.Vector3;
import org.junit.Test;
import pbd.PositionBasedDynamics;
import pbd.SPHKernel;

import org.junit.Test;
import static org.junit.Assert.*;

public class PositionBasedDynamicsTest {
	PositionBasedDynamics pbd = new PositionBasedDynamics(new SPHKernel(0));

	@Test
	public void nonStatic() {
		//PositionBasedDynamics.solveDistanceConstraint()
		//PositionBasedDynamics.solveDihedralConstraint()
		//PositionBasedDynamics.solveVolumeConstraint()
		//PositionBasedDynamics.solveEdgePointDistConstraint
		//PositionBasedDynamics.solveTrianglePointDistConstraint
		//PositionBasedDynamics.solveEdgeEdgeDistConstraint

		// iso bending
		//PositionBasedDynamics.cotTheta()
		//PositionBasedDynamics.computeQuadraticBendingMat()
		//PositionBasedDynamics.solveIsometricBendingConstraint()

		// shape matching
		//PositionBasedDynamics.computeShapeMatchingRestInfo()
		//PositionBasedDynamics.solveShapeMatchingConstraint()

		// strain based dynamics
		//PositionBasedDynamics.computeStrainTriangleInvRestMat()
		//PositionBasedDynamics.solveStrainTriangleConstraint
		//PositionBasedDynamics.computeStrainTriangleInvRestMat()
		//PositionBasedDynamics.solveStrainTetraConstraint()

		// FEM PBD
		//PositionBasedDynamics.computeFEMTetraInvRestMat()
		//PositionBasedDynamics.solveFEMTetraConstraint()
		//PositionBasedDynamics.computeFEMTriangleInvRestMat()
		//PositionBasedDynamics.solveFEMTriangleConstraint();

		// Position based fluids
		//PositionBasedDynamics.computePBFDensity()
		//PositionBasedDynamics.computePBFLagrangeMultiplier()
		//PositionBasedDynamics.solveDensityConstraint()
	}
}
