using UnityEngine;
using System.Collections;

public class CreateHypercube : MonoBehaviour {
	
	Vector4 Wa, Wb, Wc, Wd;
	Vector4 fromVec, up, over;
	float xyRot = 0f, yzRot = 0f, zxRot = 0f, xwRot = 0f, ywRot = 0f, zwRot = 0f;
	bool secondPerspective;
	MeshGenerator meshGen;
	GameObject[] leftCylinders;
	GameObject[] rightCylinders;
	
	void Start () {
		meshGen = (MeshGenerator)gameObject.GetComponent("MeshGenerator");
		leftCylinders = new GameObject[32];
		rightCylinders = new GameObject[32];
		//fromVec = new Vector4(4,0,0,-0.25f);
		fromVec = new Vector4(2.83f, 2.83f, .01f, -.03f);
		//up = new Vector4(0,1,0,0);
		up = new Vector4(-.71f, .71f, 0, 0);
		//over = new Vector4(0,0,1,0);
		over = new Vector4(0, 0, 1, .02f);
		drawHypercubeWithRotation(Matrix4x4.identity, true, leftCylinders);
		
		//fromVec = new Vector4(4,0,0,0.25f);
		fromVec = new Vector4(2.83f, 2.83f, .01f, .03f);
		drawHypercubeWithRotation(Matrix4x4.identity, false, rightCylinders);
//		for (int i = 0; i < leftCylinders.Length; i++) {
//			//leftCylinders[i].renderer.enabled = (false);
//			//rightCylinders[i].renderer.enabled = (false);
//		}
		
	}
	
	void drawHypertetrahedronWithRotation(Matrix4x4 rotationMatrix) {
		calc4Matrix();
		Vector3[,] hyperVecs = createHypertetrahedronPointsWithMatrix(rotationMatrix);
		for (int i = 0; i < hyperVecs.Length / 2; i++) {
			Vector3 point1 = hyperVecs[0,0];
			Vector3 point2 = hyperVecs[0,1];
			if (!(point1.x == point2.x && point1.y == point2.y && point1.z == point2.z)) {
				meshGen.drawCylinderWithEndpoints(hyperVecs[i,0], hyperVecs[i,1], null);	
			}
		}

	}
	
	Vector3[,] createHypertetrahedronPointsWithMatrix(Matrix4x4 rotationMatrix) {
		Vector4[,] points = new Vector4[10, 2];
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
		calc4Matrix();
		Vector3[,] hyperVecs = createHypercubePointsWithMatrix(rotationMatrix);
		for (int i = 0; i < hyperVecs.Length / 2; i++) {
			//Vector3 point1 = hyperVecs[0,0];
			//Vector3 point2 = hyperVecs[0,1];
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
			// Multiply rotation matrix to points
			point1 = multiply(matrix, point1);
			point2 = multiply(matrix, point2);
			hyperVecs[k,0] = projectTo3D(point1);
			hyperVecs[k,1] = projectTo3D(point2);
		}
		return hyperVecs;
	}
	
	string[,] removeDuplicates(string[,] strings) {
		// REMOVE DUPLICATES
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
		float vAngle = Mathf.PI / 4.0f;
		// Perspective projection
		t = 1.0f / Mathf.Tan(vAngle / 2.0f);
		
		v = inputVec - fromVec;
		s = t / Vector4.Dot(v, Wd);
		projectedPoint.x = s * Vector4.Dot(v, Wa);
		projectedPoint.y = s * Vector4.Dot(v, Wb);
		projectedPoint.z = s * Vector4.Dot(v, Wc);
		return projectedPoint;
	}
	
	void calc4Matrix() {
		// Calculates the 4d perspective matrix to apply to the hyperpoints
		Vector4 toVec = new Vector4(0,0,0,0);
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
		Vector4 product = new Vector4();
		product.x = a.y*(b.z*c.w - c.z*b.w) - a.z*(b.y*c.w - c.y*b.w) + a.w*(b.y*c.z - c.y*b.z);
		product.y = -a.x*(b.z*c.w - c.z*b.w) + a.z*(b.x*c.w - c.x*b.w) - a.w*(b.x*c.z - c.x*b.z);
		product.z = a.x*(b.y*c.w - c.y*b.w) - a.y*(b.x*c.w - c.x*b.w) + a.w*(b.x*c.y - c.x*b.y);
		product.w = -a.x*(b.y*c.z - c.y*b.z) + a.y*(b.x*c.z - c.x*b.z) - a.z*(b.x*c.y - c.x*b.y);
		return product;
	}
		
	void FixedUpdate() {
//		xyRot +=  System.Convert.ToInt32(Input.GetKey(KeyCode.Q)) * Mathf.Deg2Rad;
//		Matrix4x4 xyRotMatrix = xyRotationBy(xyRot);
//		
//		yzRot += System.Convert.ToInt32(Input.GetKey(KeyCode.W)) * Mathf.Deg2Rad;
//		Matrix4x4 yzRotMatrix = yzRotationBy(yzRot);
//		
//		zxRot += System.Convert.ToInt32(Input.GetKey(KeyCode.E)) * Mathf.Deg2Rad;
//		Matrix4x4 zxRotMatrix = zxRotationBy(zxRot);
//		
//		xwRot += System.Convert.ToInt32(Input.GetKey(KeyCode.R)) * Mathf.Deg2Rad;
//		Matrix4x4 xwRotMatrix = xwRotationBy(xwRot);
//		
//		ywRot += System.Convert.ToInt32(Input.GetKey(KeyCode.T)) * Mathf.Deg2Rad;
//		Matrix4x4 ywRotMatrix = ywRotationBy(ywRot);
//		
//		zwRot += System.Convert.ToInt32(Input.GetKey(KeyCode.Y)) * Mathf.Deg2Rad;
//		Matrix4x4 zwRotMatrix = zwRotationBy(zwRot);
//
//		Matrix4x4 combinedMatrix = multiply(multiply(multiply(zxRotMatrix, xwRotMatrix), multiply(xyRotMatrix, yzRotMatrix)), multiply(ywRotMatrix, zwRotMatrix));
		//drawHypertetrahedronWithRotation(combinedMatrix);
		//drawHypercubeWithRotation(combinedMatrix);
		//meshGen.nextIndex = 0;
		
	}
	
	public void updateRotationsWithMove(Quaternion move) {
		Matrix4x4 rotation = Matrix4x4.TRS(Vector3.zero, move, Vector3.one);
		fromVec = multiply(rotation, new Vector4(2.83f, 2.83f, .01f, -.03f));
		Debug.Log(fromVec);
		drawHypercubeWithRotation(Matrix4x4.identity, true, leftCylinders);
		fromVec = multiply(rotation, new Vector4(2.83f, 2.83f, .01f, .03f));
		drawHypercubeWithRotation(Matrix4x4.identity, false, rightCylinders);

	}
	public void OnWillRenderObject(){
		if (Camera.current.name == "CameraLeft") {
			for (int i = 0; i < rightCylinders.Length; i++) {
				leftCylinders[i].renderer.material.SetFloat("_Show",0.0f);
				rightCylinders[i].renderer.material.SetFloat("_Show",1.0f);
			}
		} else if (Camera.current.name == "CameraRight") {
			for (int i = 0; i < rightCylinders.Length; i++) {
				leftCylinders[i].renderer.material.SetFloat("_Show",1.0f);
				rightCylinders[i].renderer.material.SetFloat("_Show",0.0f);
			}
		}	
	}
	
	void Update () {
		if (Input.GetKeyDown(KeyCode.P)) {
			togglePerspective();	
		}
	}
	
	void togglePerspective() {
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
		Vector4 returnVal = new Vector4();
		Vector4 row0 = mat.GetRow(0); Vector4 row1 = mat.GetRow(1); Vector4 row2 = mat.GetRow(2); Vector4 row3 = mat.GetRow(3);
		returnVal.x = row0.x*point.x + row0.y*point.x + row0.z*point.x + row0.w*point.x;
		returnVal.y = row1.x*point.y + row1.y*point.y + row1.z*point.y + row1.w*point.y;
		returnVal.z = row2.x*point.z + row2.y*point.z + row2.z*point.z + row2.w*point.z;
		returnVal.w = row3.x*point.w + row3.y*point.w + row3.z*point.w + row3.w*point.w;
		return returnVal;
	}
}
