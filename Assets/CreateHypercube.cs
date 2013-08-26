using UnityEngine;
using System.Collections;

public class CreateHypercube : MonoBehaviour {
	
	Vector4 Wa, Wb, Wc, Wd;
	Vector4 fromVec, up, over;
	Vector4 toVec;
	float xyRot = 0f, yzRot = 0f, zxRot = 0f, xwRot = 0f, ywRot = 0f, zwRot = 0f;
	bool secondPerspective;
	MeshGenerator meshGen;
	GameObject[] leftCylinders;
	GameObject[] rightCylinders;
	public bool usingOculus;
	public bool hypercube;
	private Quaternion oculusMove;
	
	void Start () {
		meshGen = (MeshGenerator)gameObject.GetComponent("MeshGenerator");
		leftCylinders = new GameObject[(hypercube) ? 32 : 10];
		rightCylinders = new GameObject[(hypercube) ? 32 : 10];
		
		// Projection values pulled from http://steve.hollasch.net/thesis/chapter4.html
		toVec = new Vector4(1,0,0,0);
		fromVec = new Vector4(0,0,0,-.03f);
		up = new Vector4(0,1,0,0);
		over = new Vector4(0,0,1,0);
		if (hypercube) {
			drawHypercubeWithRotation(Matrix4x4.identity, true, leftCylinders);
		} else {
			drawHypertetrahedronWithRotation(Matrix4x4.identity, true, leftCylinders);
		}
		
		if (usingOculus) {
			fromVec = new Vector4(0,0,0,.03f);
			if (hypercube) {
				drawHypercubeWithRotation(Matrix4x4.identity, false, rightCylinders);
			} else {
				drawHypertetrahedronWithRotation(Matrix4x4.identity, false, rightCylinders);
			}
		}
	}
	
	void drawHypertetrahedronWithRotation(Matrix4x4 rotationMatrix, bool leftProjection, GameObject[] gameObjects) {
		calc4Matrix(leftProjection);
		Vector3[,] hyperVecs = createHypertetrahedronPointsWithMatrix(rotationMatrix);
		for (int i = 0; i < hyperVecs.Length / 2; i++) {
			GameObject go = meshGen.drawCylinderWithEndpoints(hyperVecs[i,0], hyperVecs[i,1], gameObjects[i]);
			if (leftProjection) {
				leftCylinders[i] = go;	
			} else {
				rightCylinders[i] = go;
			}
		}

	}
	
	Vector3[,] createHypertetrahedronPointsWithMatrix(Matrix4x4 rotationMatrix) {
		Vector4[,] points = new Vector4[10, 2];
		// Vertices obtained at http://en.wikipedia.org/wiki/5-cell
		points[0,0] = new Vector4(1,1,1,-1); points[0,1] = new Vector4(1,-1,-1,-1);
		points[1,0] = new Vector4(1,-1,-1,-1); points[1,1] = new Vector4(-1,1,-1,-1);
		points[2,0] = new Vector4(-1,1,-1,-1); points[2,1] = new Vector4(-1,-1,1,-1);
		points[3,0] = new Vector4(1,1,1,-1); points[3,1] = new Vector4(-1,-1,1,-1);
		points[4,0] = new Vector4(1,1,1,-1); points[4,1] = new Vector4(-1,1,-1,-1);
		points[5,0] = new Vector4(-1,-1,1,-1); points[5,1] = new Vector4(1,-1,-1,-1);
		points[6,0] = new Vector4(1,1,1,-1); points[6,1] = new Vector4(0,0,0,Mathf.Sqrt(5) - 1);
		points[7,0] = new Vector4(1,-1,-1,-1); points[7,1] = new Vector4(0,0,0,Mathf.Sqrt(5) - 1);
		points[8,0] = new Vector4(-1,1,-1,-1); points[8,1] = new Vector4(0,0,0,Mathf.Sqrt(5) - 1);
		points[9,0] = new Vector4(-1,-1,1,-1); points[9,1] = new Vector4(0,0,0,Mathf.Sqrt(5) - 1);
		Vector3[,]hyperVecs = new Vector3[points.Length/2, 2];
		for (int i = 0; i < points.Length/2; i++) {
			Vector4 point1 = points[i,0];
			Vector4 point2 = points[i,1];
			point1 = multiply(rotationMatrix, point1);
			point2 = multiply(rotationMatrix, point2);
			hyperVecs[i,0] = projectTo3D(point1);
			hyperVecs[i,1] = projectTo3D(point2);
		}
		return hyperVecs;
	}
	void drawHypercubeWithRotation(Matrix4x4 rotationMatrix, bool leftProjection, GameObject[] gameObjects) {
		calc4Matrix(leftProjection);
		Vector3[,] hyperVecs = createHypercubePointsWithMatrix(rotationMatrix);
		for (int i = 0; i < hyperVecs.Length / 2; i++) {
			GameObject go = meshGen.drawCylinderWithEndpoints(hyperVecs[i,0], hyperVecs[i,1], gameObjects[i]);	
			if (leftProjection) {
				leftCylinders[i] = go;	
			} else {
				rightCylinders[i] = go;
			}
		}
	}
	string binaryValue(int input) {
		char[] binArray = System.Convert.ToString(input, 2).ToCharArray();
		if (binArray.Length < 4) {
			string padding = "";
			for (int i = 0; i < 4 - binArray.Length; i++) {
				padding	= System.String.Concat("0", padding);
			}
				
			return System.String.Concat(padding, new System.String(binArray));
		}
		return new System.String(binArray);
	
	}
	
	int intValue(string binary) {
		return System.Convert.ToInt32(binary, 2);	
	}
	
	string changeBinaryAtIndex(string binVal, int idx) {
		char[] chars = binVal.ToCharArray();
		if (chars[idx] == '0') {
			chars[idx] = '1';	
		} else {
			chars[idx] = '0';	
		}
		return new System.String(chars);
	}
	
	Vector3[,] createHypercubePointsWithMatrix(Matrix4x4 matrix) {
		string[,] strings = new string[64,2];
		for (int i = 0; i <= 15; i++) {
			string binVal = binaryValue(i);
			for (int k = 0; k < 4; k++) {
				string newVal = changeBinaryAtIndex(binVal, k);	
				strings[i * 4 + k,0] = binVal;
				strings[i * 4 + k,1] = newVal;

			}
		}
		string[,] hyperPoints = removeDuplicates(strings);
		Vector3[,] hyperVecs = new Vector3[hyperPoints.Length/2,2];
		for (int k = 0; k < hyperPoints.Length / 2; k++) {
			Vector4 point1 = convertToPoint(hyperPoints[k,0]);
			Vector4 point2 = convertToPoint(hyperPoints[k,1]);
			// Transform Points to (+-1,+-1,+-1,+-1);
			point1 = point1 * 2.0f; point1.x -= 1; point1.y -= 1; point1.z -=1; point1.w -=1;
			point2 = point2 * 2.0f; point2.x -=1; point2.y -=1; point2.z -=1; point2.w -=1;
			point1 = multiply(matrix, point1);
			point2 = multiply(matrix, point2);
			hyperVecs[k,0] = projectTo3D(point1);
			hyperVecs[k,1] = projectTo3D(point2);
		}
		return hyperVecs;
	}
	
	string[,] removeDuplicates(string[,] strings) {
		string[,] hyperPoints = new string[32,2];
		int nextIndex = 0;
		for (int j = 0; j < strings.Length / 2; j++) {
			bool newVal = true;
			for (int l = 0; l < hyperPoints.Length / 2; l++) {
				if (j != l &&(hyperPoints[l,0] != "" && hyperPoints[l,1] != "")) {
					string compString1 = System.String.Concat(strings[j,0], strings[j,1]);
					string compString2 = System.String.Concat(hyperPoints[l,0], hyperPoints[l,1]);
					string reverseCompString2 = System.String.Concat(hyperPoints[l,1], hyperPoints[l,0]);
					if (System.String.Compare(compString1, compString2, true) == 0 || System.String.Compare(compString1, reverseCompString2) == 0) {
						newVal = false;
					}
				}
			}
			if (newVal) {
				hyperPoints[nextIndex,0] = strings[j,0];
				hyperPoints[nextIndex++,1] = strings[j,1];
			}
		}	
		return hyperPoints;
	}
	
	Vector4 convertToPoint(string inputString) {
		Vector4 vec = new Vector4(System.Convert.ToInt32(inputString[0].ToString()), System.Convert.ToInt32(inputString[1].ToString()), System.Convert.ToInt32(inputString[2].ToString()), System.Convert.ToInt32(inputString[3].ToString()));
		return vec;
	}
	
	Vector3 projectTo3D(Vector4 inputVec) {
		Vector3 projectedPoint = new Vector3();
		float s,t;
		Vector4 v;
		float vAngle = Mathf.PI / 10.0f;
		// Perspective projection
		t = 1.0f / Mathf.Tan(vAngle / 2.0f);
		
		v = inputVec - fromVec;
		s = t / Vector4.Dot(v, Wd);
		projectedPoint.x = s * Vector4.Dot(v, Wa);
		projectedPoint.y = s * Vector4.Dot(v, Wb);
		projectedPoint.z = s * Vector4.Dot(v, Wc);
		return projectedPoint;
	}
	
	void calc4Matrix(bool leftProjection) {
		// Calculates the 4d perspective matrix to apply to the 4D points
		float norm;
		
	 	// Get the normalized Wd column-vector.
		Wd = toVec - fromVec;
		norm = Wd.magnitude;
		Wd.Scale(new Vector4(1.0f/norm,1.0f/norm,1.0f/norm,1.0f/norm));
		
		// Calculate the normalized Wa column-vector.
		Wa = cross(up, over, Wd);
		norm = Wa.magnitude;
		Wa.Scale(new Vector4(1.0f/norm,1.0f/norm,1.0f/norm,1.0f/norm));
		
		//Calculate the normalized Wb column-vector
		Wb = cross(over, Wd, Wa);
		norm = Wb.magnitude;
		Wb.Scale(new Vector4(1.0f/norm,1.0f/norm,1.0f/norm,1.0f/norm));
		
		// Calculate the Wc column-vector.
		Wc = cross(Wd, Wa, Wb);
	}
	
	Vector4 cross(Vector4 a, Vector4 b, Vector4 c) {
		// Calculates the cross product of three Vector4's
		Vector4 product = new Vector4();
		product.x = a.y*(b.z*c.w - c.z*b.w) - a.z*(b.y*c.w - c.y*b.w) + a.w*(b.y*c.z - c.y*b.z);
		product.y = -a.x*(b.z*c.w - c.z*b.w) + a.z*(b.x*c.w - c.x*b.w) - a.w*(b.x*c.z - c.x*b.z);
		product.z = a.x*(b.y*c.w - c.y*b.w) - a.y*(b.x*c.w - c.x*b.w) + a.w*(b.x*c.y - c.x*b.y);
		product.w = -a.x*(b.y*c.z - c.y*b.z) + a.y*(b.x*c.z - c.x*b.z) - a.z*(b.x*c.y - c.x*b.y);
		return product;
	}
		
	void FixedUpdate() {
		if (!usingOculus) {
			// Calculate total 4d rotation from keys 'Q,W,E,R,T,Y' mapping to 'xy,yz,xz,xw,yw,zw' rotations
			xyRot +=  System.Convert.ToInt32(Input.GetKey(KeyCode.Q)) * Mathf.Deg2Rad;
			Matrix4x4 xyRotMatrix = xyRotationBy(xyRot);
			
			yzRot += System.Convert.ToInt32(Input.GetKey(KeyCode.W)) * Mathf.Deg2Rad;
			Matrix4x4 yzRotMatrix = yzRotationBy(yzRot);
			
			zxRot += System.Convert.ToInt32(Input.GetKey(KeyCode.E)) * Mathf.Deg2Rad;
			Matrix4x4 zxRotMatrix = zxRotationBy(zxRot);
			
			xwRot += System.Convert.ToInt32(Input.GetKey(KeyCode.R)) * Mathf.Deg2Rad;
			Matrix4x4 xwRotMatrix = xwRotationBy(xwRot);
			
			ywRot += System.Convert.ToInt32(Input.GetKey(KeyCode.T)) * Mathf.Deg2Rad;
			Matrix4x4 ywRotMatrix = ywRotationBy(ywRot);
			
			zwRot += System.Convert.ToInt32(Input.GetKey(KeyCode.Y)) * Mathf.Deg2Rad;
			Matrix4x4 zwRotMatrix = zwRotationBy(zwRot);
			Matrix4x4 combinedMatrix = multiply(multiply(multiply(zxRotMatrix, xwRotMatrix), multiply(xyRotMatrix, yzRotMatrix)), multiply(ywRotMatrix, zwRotMatrix));
			if (hypercube) {
				drawHypercubeWithRotation(combinedMatrix, true, leftCylinders);
			} else {
				drawHypertetrahedronWithRotation(combinedMatrix, true, leftCylinders);
			}
		}		
	}
	
	public void updateRotationsWithMove(Quaternion move) { 
		// This method is utilized only for the Oculus integration

		// Ensure a valid Quaternion value
		oculusMove = move;
		float sum = 0;
    	for (int i = 0; i < 4; ++i) {
        	sum += oculusMove[i] * oculusMove[i];
		}
		float magnitudeInverse = 1 / Mathf.Sqrt(sum);
		for (int i = 0; i < 4; ++i) {
			 oculusMove[i] *= magnitudeInverse;  
		}
		
		float xAxis = halveAngle(oculusMove.eulerAngles.x);
		float yAxis = halveAngle(oculusMove.eulerAngles.y);
		float zAxis = halveAngle(oculusMove.eulerAngles.z);
				
		Matrix4x4 rotation = xyRotationBy(Mathf.Deg2Rad * zAxis);
		rotation = multiply(rotation, zxRotationBy(Mathf.Deg2Rad * yAxis));
		rotation = multiply(rotation, yzRotationBy(- Mathf.Deg2Rad * xAxis));
		
		Vector4 moddedFromVec = multiply(rotation, new Vector4(-.03f, 0, 4, 0));
		fromVec = reorderCoordinates(moddedFromVec);
		Vector4 moddedUp = multiply(rotation, new Vector4(0,1,0,0));
		up = reorderCoordinates(moddedUp);
		Vector4 moddedOver = multiply(rotation, new Vector4(0,0,0,1));
		over = reorderCoordinates(moddedOver);
		
		Matrix4x4 trs = Matrix4x4.TRS(Vector3.zero, Quaternion.identity, new Vector3(1,1,1));
		if (hypercube) {
			drawHypercubeWithRotation(trs, true, leftCylinders);
		} else {
			drawHypertetrahedronWithRotation(Matrix4x4.identity, true, leftCylinders);
		}
		moddedFromVec = multiply(rotation, new Vector4(.03f, 0, 4, 0));
		fromVec = reorderCoordinates(moddedFromVec);

		if (hypercube) {
			drawHypercubeWithRotation(trs, false, rightCylinders);
		} else {
			drawHypertetrahedronWithRotation(Matrix4x4.identity, false, rightCylinders);
		}

	}
	
	private float halveAngle(float ang) {
		if (ang >= 180.0f) {
			return (ang - 360.0f) * 0.5f;	
		} else {
			return ang * 0.5f;
		}
	}
	
	private Vector4 reorderCoordinates(Vector4 inp) {
		Vector4 reorderedVec = new Vector4();
		reorderedVec.w = inp.x;
		reorderedVec.y = inp.y;
		reorderedVec.z = inp.w;
		reorderedVec.x = inp.z;
		return reorderedVec;
	}
	
	public void restrictFromVector() {
		// Restricts the fromVector so that the hyperObject doesn't flip or stretch abnormally
		fromVec.x = Mathf.Abs(fromVec.x);
		fromVec.y = Mathf.Abs(fromVec.y);
		if (fromVec.x + fromVec.y < 3.5f) {
			if (fromVec.x < fromVec.y) {
				fromVec.x = 3.5f - fromVec.y;	
			} else {
				fromVec.y = 3.5f - fromVec.x;	
			}
		}
	}
	public void OnWillRenderObject(){ 
		// This method is utilized only for the Oculus integration
		// Toggles the different 4d projections based on which camera is being displayed.
		if (Camera.current.name == "CameraLeft") {
			for (int i = 0; i < rightCylinders.Length; i++) {
				leftCylinders[i].renderer.material.SetFloat("_Show",1.0f);
				rightCylinders[i].renderer.material.SetFloat("_Show",0.0f);
			}
		} else if (Camera.current.name == "CameraRight") {
			for (int i = 0; i < rightCylinders.Length; i++) {
				leftCylinders[i].renderer.material.SetFloat("_Show",0.0f);
				rightCylinders[i].renderer.material.SetFloat("_Show",1.0f);
			}
		}	
	}
	
	void Update () {
		if (Input.GetKeyDown(KeyCode.P)) {
			togglePerspective();	
		}
	}
	
	void togglePerspective() { 
		// Toggles perspective from keyboard input 'P'. 
		// Switches between two perspectives found at: http://steve.hollasch.net/thesis/chapter4.html
		secondPerspective = !secondPerspective;
		if (secondPerspective) {
			fromVec = new Vector4(2.83f, 2.83f, .01f, 0);
			up = new Vector4(-.71f, .71f, 0, 0);
			over = new Vector4(0, 0, 1, .02f);
		} else {
			fromVec = new Vector4(4,0,0,0);
			up = new Vector4(0,1,0,0);
			over = new Vector4(0,0,1,0);
		}
	
	}
	
	// 4D Rotation utility methods
	Matrix4x4 xyRotationBy(float radians) {
		Matrix4x4 rotMatrix = new Matrix4x4();
		rotMatrix.SetRow(0, new Vector4(Mathf.Cos(radians), Mathf.Sin(radians), 0, 0));
		rotMatrix.SetRow(1, new Vector4(-Mathf.Sin(radians), Mathf.Cos(radians), 0, 0));
		rotMatrix.SetRow(2, new Vector4(0, 0, 1, 0));
		rotMatrix.SetRow(3, new Vector4(0, 0, 0, 1));
		return rotMatrix;
	}
	
	Matrix4x4 yzRotationBy(float radians) {
		Matrix4x4 rotMatrix = new Matrix4x4();
		rotMatrix.SetRow(0, new Vector4(1, 0, 0, 0));
		rotMatrix.SetRow(1, new Vector4(0, Mathf.Cos(radians), Mathf.Sin(radians), 0));
		rotMatrix.SetRow(2, new Vector4(0, -Mathf.Sin(radians), Mathf.Cos(radians), 0));
		rotMatrix.SetRow(3, new Vector4(0, 0, 0, 1));
		return rotMatrix;
	}
	
	Matrix4x4 zxRotationBy(float radians) {
		Matrix4x4 rotMatrix = new Matrix4x4();
		rotMatrix.SetRow(0, new Vector4(Mathf.Cos(radians), 0, -Mathf.Sin(radians), 0));
		rotMatrix.SetRow(1, new Vector4(0, 1, 0, 0));
		rotMatrix.SetRow(2, new Vector4(Mathf.Sin(radians), 0, Mathf.Cos(radians), 0));
		rotMatrix.SetRow(3, new Vector4(0, 0, 0, 1));
		return rotMatrix;	
	}
	
	Matrix4x4 xwRotationBy(float radians) {
		Matrix4x4 rotMatrix = new Matrix4x4();
		rotMatrix.SetRow(0, new Vector4(Mathf.Cos(radians), 0, Mathf.Sin(radians), 0));
		rotMatrix.SetRow(1, new Vector4(0, 1, 0, 0));
		rotMatrix.SetRow(2, new Vector4(0, 0, 1, 0));
		rotMatrix.SetRow(3, new Vector4(-Mathf.Sin(radians), 0, 0, Mathf.Cos(radians)));
		return rotMatrix;	
	}
	
	Matrix4x4 ywRotationBy(float radians) {
		Matrix4x4 rotMatrix = new Matrix4x4();
		rotMatrix.SetRow(0, new Vector4(1, 0, 0, 0));
		rotMatrix.SetRow(1, new Vector4(0, Mathf.Cos(radians), 0, -Mathf.Sin(radians)));
		rotMatrix.SetRow(2, new Vector4(0, 0, 1, 0));
		rotMatrix.SetRow(3, new Vector4(0, Mathf.Sin(radians), 0, Mathf.Cos(radians)));
		return rotMatrix;	
	}
	
	Matrix4x4 zwRotationBy(float radians) {
		Matrix4x4 rotMatrix = new Matrix4x4();
		rotMatrix.SetRow(0, new Vector4(1, 0, 0, 0));
		rotMatrix.SetRow(1, new Vector4(0, 1, 0, 0));
		rotMatrix.SetRow(2, new Vector4(0, 0, Mathf.Cos(radians), -Mathf.Sin(radians)));
		rotMatrix.SetRow(3, new Vector4(0, 0, Mathf.Sin(radians), Mathf.Cos(radians)));
		return rotMatrix;	
	}
	
	Matrix4x4 multiply(Matrix4x4 src1, Matrix4x4 src2) {
		// Multiplies two 4x4 matrices
		Matrix4x4 dest = new Matrix4x4();
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				dest[i,j] = 0;
				for (int k = 0; k < 4; k++) {
					dest[i,j] += src1[i,k] * src2[k,j];
				}
			}
		}
		return dest;
	}
	
	Vector4 multiply(Matrix4x4 mat, Vector4 point) {
		// Multiplies a 4x4 matrix by a vector4 
		Vector4 returnVal = new Vector4();
		Vector4 row0 = mat.GetRow(0); Vector4 row1 = mat.GetRow(1); Vector4 row2 = mat.GetRow(2); Vector4 row3 = mat.GetRow(3);
		returnVal.x = row0.x*point.x + row0.y*point.y + row0.z*point.z + row0.w*point.w;
		returnVal.y = row1.x*point.x + row1.y*point.y + row1.z*point.z + row1.w*point.w;
		returnVal.z = row2.x*point.x + row2.y*point.y + row2.z*point.z + row2.w*point.w;
		returnVal.w = row3.x*point.x + row3.y*point.y + row3.z*point.z + row3.w*point.w;
		return returnVal;
	}
	
	Vector4 multiply(Vector4 point, Matrix4x4 mat) {
		Vector4 returnVal = new Vector4();
		Vector4 col0 = mat.GetColumn(0); Vector4 col1 = mat.GetColumn(1); Vector4 col2 = mat.GetColumn(2); Vector4 col3 = mat.GetColumn(3);
		returnVal.x = col0.x*point.x + col0.y*point.y + col0.z*point.z + col0.w*point.w;
		returnVal.y = col1.x*point.x + col1.y*point.y + col1.z*point.z + col1.w*point.w;
		returnVal.z = col2.x*point.x + col2.y*point.y + col2.z*point.z + col2.w*point.w;
		returnVal.w = col3.x*point.x + col3.y*point.y + col3.z*point.z + col3.w*point.w;
		return returnVal;
	}
}
