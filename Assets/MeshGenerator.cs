using UnityEngine;
using System.Collections;

public class MeshGenerator : MonoBehaviour {
	public GameObject cylinderPrefab;
	Vector3[] vertices;
	int[] triangles;
	Vector3[] normals;
	int trianglesIndex = 0;
	Vector2[] uvs;
	public Transform target;
	GameObject[] cylinders;
	public int nextIndex = 0;

	// Use this for initialization
	void Start () {		
		cylinders = new GameObject[50];
		//drawCylinderWithEndpoints(new Vector3(0,0,0), new Vector3(1,1,1));
	}
	
	public void drawCylinderWithEndpoints(Vector3 startPoint, Vector3 endPoint) {
		if (startPoint != endPoint) {
			float xD = endPoint.x - startPoint.x;
			float yD = endPoint.y - startPoint.y;
			float zD = endPoint.z - startPoint.z;
			float distance = Mathf.Sqrt(xD*xD + yD*yD + zD*zD);
			
			if (GameObject.Find(System.String.Format("{0}", nextIndex))) {
				cylinderPrefab = GameObject.Find(System.String.Format("{0}", nextIndex));
				cylinderPrefab.transform.position = Vector3.zero;
				cylinderPrefab.transform.rotation = Quaternion.identity;
			} else {
				cylinderPrefab = (GameObject)Instantiate(cylinderPrefab, Vector3.zero, Quaternion.identity);
			}
			createCylinder(0.05f, distance, 8, cylinderPrefab);
			cylinderPrefab.transform.Translate(startPoint);
			cylinderPrefab.transform.LookAt(endPoint);	
			cylinderPrefab.name = System.String.Format("{0}", nextIndex);
			cylinders[nextIndex++] = cylinderPrefab;
		}			 
	}
	
	public void destroyCylinders() {
		for (int i = 0; i < cylinders.Length; i++) {
			GameObject.Destroy(cylinders[i]);
		}
	}
	
	float calculateRotation(Vector2 startPoint, Vector2 endPoint) {
		Vector2 legStartPoint = new Vector2(endPoint.x, startPoint.y);
		Vector2 legEndPoint = new Vector2(endPoint.x, endPoint.y);
		float legDistance = Mathf.Sqrt(Mathf.Pow((legEndPoint.x - legStartPoint.x), 2) + Mathf.Pow((legEndPoint.y - legStartPoint.y), 2));
		float hypDistance = Mathf.Sqrt(Mathf.Pow((endPoint.x - startPoint.x), 2) + Mathf.Pow((endPoint.y - startPoint.y), 2));
		
		float angle = Mathf.Asin(legDistance / hypDistance) * (180.0f / Mathf.PI);
		return angle;
	}
	
	void createWireframe() {
		Vector3[] rotTrans1 = new Vector3[]{new Vector3(0,0,0), new Vector3(0,0,0)};
		Vector3[] rotTrans2 = new Vector3[]{new Vector3(1,0,0), new Vector3(0,0,0)};
		Vector3[] rotTrans3 = new Vector3[]{new Vector3(1,0,0), new Vector3(0,-1,0)};
		Vector3[] rotTrans4 = new Vector3[]{new Vector3(0,0,0), new Vector3(0,0,-1)};
		
		Vector3[] rotTrans5 = new Vector3[]{rotTrans1[0], new Vector3(1,0,0)};
		Vector3[] rotTrans6 = new Vector3[]{rotTrans2[0], new Vector3(1,0,0)};
		Vector3[] rotTrans7 = new Vector3[]{rotTrans3[0], new Vector3(1,-1,0)};
		Vector3[] rotTrans8 = new Vector3[]{rotTrans4[0], new Vector3(1,0,-1)};
		
		Vector3[] rotTrans9 = new Vector3[]{new Vector3(0,0,1), new Vector3(0,-1,0)};
		Vector3[] rotTrans10 = new Vector3[]{rotTrans9[0], new Vector3(0,-1,-1)};
		Vector3[] rotTrans11 = new Vector3[]{rotTrans9[0], new Vector3(0,0,0)};
		Vector3[] rotTrans12 = new Vector3[]{rotTrans9[0], new Vector3(0,0,-1)};
		
		Vector3[] rotTrans13 = new Vector3[]{rotTrans2[0], rotTrans2[1]};
		Vector3[] rotTrans14 = new Vector3[]{rotTrans6[0], rotTrans6[1]};
		Vector3[] rotTrans15 = new Vector3[]{rotTrans10[0], rotTrans10[1]};
		Vector3[] rotTrans16 = new Vector3[]{new Vector3(0,0,1),  rotTrans6[1]};
		Vector3[] rotTrans17 = new Vector3[]{new Vector3(0,0,1), rotTrans2[1]};
		Vector3[] rotTrans18 = new Vector3[]{new Vector3(0,0,1), new Vector3(0,-1,0)};
		
		Vector3[][] rotTrans = new Vector3[][]{rotTrans1, rotTrans2, rotTrans3, rotTrans4,
											   rotTrans5, rotTrans6, rotTrans7, rotTrans8,
											   rotTrans9, rotTrans10, rotTrans11, rotTrans12,
											   rotTrans13, rotTrans14, rotTrans15, rotTrans16, rotTrans17, rotTrans18};
		
		for (int i = 0; i < 18; i++) {
			int cylName = i + 1;
			float rot = 90.0f;
			Matrix4x4 mat = Matrix4x4.identity;
			if (i >= 12) {
				if (i < 14 || i > 15) {
					rot = 45.0f;
				} else if (i < 15) {
					rot = 135f;
				} else {
					rot = -45f;	
				}
				mat = Matrix4x4.Scale(new Vector3(.1f,Mathf.Sqrt(2.0f),1.0f));
			}
			GameObject cyl = GameObject.Find(string.Format("Cylinder{0}", cylName));
			createCylinder(0.1f,1.0f,8,cyl,mat);
			cyl.transform.Translate(rotTrans[i][1]);

			if (i >15) {
				cyl.transform.Rotate(new Vector3(1,0,0), 90.0f);
			}
			cyl.transform.Rotate(rotTrans[i][0], rot);
		}	
		
		//Scaled X's
	}

	void createTriangle() {
		Mesh triangleMesh = new Mesh();
		Vector3[] triangleVerts = new Vector3[3];
		triangleVerts[0] = new Vector3(-1,0,0);
		triangleVerts[1] = new Vector3(1,0,0);
		triangleVerts[2] = new Vector3(0,1,0);
		triangleMesh.vertices = triangleVerts;
		
		int[] triangles = new int[]{0,2,1};
		triangleMesh.triangles = triangles;
		
		Vector3[] norms = new Vector3[3];
		norms[0] = Vector3.Normalize(Vector3.Cross(triangleVerts[2] - triangleVerts[1], triangleVerts[2] - triangleVerts[0]));
		norms[1] = Vector3.Normalize(Vector3.Cross(triangleVerts[0] - triangleVerts[1], triangleVerts[2] - triangleVerts[1]));
		norms[2] = Vector3.Normalize(Vector3.Cross(triangleVerts[1] - triangleVerts[2], triangleVerts[0] - triangleVerts[2]));
		triangleMesh.normals = norms;
		
		Vector2[] uvs = new Vector2[3];
		uvs[0] = new Vector2(0,0);
		uvs[1] = new Vector2(1,0);
		uvs[2] = new Vector2(0.5f,1.0f);
		triangleMesh.uv = uvs;
		
		MeshFilter filter = (MeshFilter)gameObject.GetComponent("MeshFilter");
		filter.mesh = triangleMesh;
	}
	
	void createTetrahedron() {
		Mesh tetraMesh = new Mesh();
		Vector3[] tetraVerts = new Vector3[12];
		tetraVerts[0] = new Vector3(-1,0,0);
		tetraVerts[1] = new Vector3(1,0,0);
		tetraVerts[2] = new Vector3(0,1,0);
		
		tetraVerts[3] = new Vector3(0,0,1);
		tetraVerts[4] = new Vector3(-1,0,0);
		tetraVerts[5] = new Vector3(0,1,0);
		
		tetraVerts[6] = new Vector3(0,0,1);
		tetraVerts[7] = new Vector3(0,1,0);
		tetraVerts[8] = new Vector3(1,0,0);
		
		tetraVerts[9] = new Vector3(-1,0,0);
		tetraVerts[10] = new Vector3(1,0,0);
		tetraVerts[11] = new Vector3(0,0,1);
		tetraMesh.vertices = tetraVerts;
		
		int[] triangles = new int[]{0,2,1,3,5,4,8,7,6,9,10,11};
		tetraMesh.triangles = triangles;
		
		Vector3[] norms = new Vector3[12];
		norms[0] = Vector3.Normalize(Vector3.Cross(tetraVerts[2] - tetraVerts[0], tetraVerts[1] - tetraVerts[0]));
		norms[1] = Vector3.Normalize(Vector3.Cross(tetraVerts[0] - tetraVerts[1], tetraVerts[2] - tetraVerts[1]));
		norms[2] = Vector3.Normalize(Vector3.Cross(tetraVerts[1] - tetraVerts[2], tetraVerts[0] - tetraVerts[2]));
		
		norms[3] = Vector3.Normalize(Vector3.Cross(tetraVerts[4] - tetraVerts[3], tetraVerts[5] - tetraVerts[3]));
		norms[4] = Vector3.Normalize(Vector3.Cross(tetraVerts[5] - tetraVerts[4], tetraVerts[3] - tetraVerts[4]));
		norms[5] = Vector3.Normalize(Vector3.Cross(tetraVerts[3] - tetraVerts[5], tetraVerts[4] - tetraVerts[5]));
		
		norms[6] = Vector3.Normalize(Vector3.Cross(tetraVerts[7] - tetraVerts[6], tetraVerts[8] - tetraVerts[6]));
		norms[7] = Vector3.Normalize(Vector3.Cross(tetraVerts[8] - tetraVerts[7], tetraVerts[6] - tetraVerts[7]));
		norms[8] = Vector3.Normalize(Vector3.Cross(tetraVerts[6] - tetraVerts[8], tetraVerts[7] - tetraVerts[8]));
		
		norms[9] = Vector3.Normalize(Vector3.Cross(tetraVerts[11] - tetraVerts[9], tetraVerts[10] - tetraVerts[9]));
		norms[10] = Vector3.Normalize(Vector3.Cross(tetraVerts[9] - tetraVerts[10], tetraVerts[11] - tetraVerts[10]));
		norms[11] = Vector3.Normalize(Vector3.Cross(tetraVerts[10] - tetraVerts[11], tetraVerts[9] - tetraVerts[11]));
			
		tetraMesh.normals = norms;
		
		Vector2[] uvs = new Vector2[12];
		Vector2 left = new Vector2(0.0f,0.0f);
		Vector2 right = new Vector2(1.0f,0.0f);
		Vector2 top = new Vector2(0.5f,1.0f);
		
		for( int i = 0; i < 12; i+=3){
			uvs[i+0] = new Vector2(left.x,left.y);
			uvs[i+1] = new Vector2(left.x,left.y);
			uvs[i+2] = new Vector2(top.x,top.y);
		}
		
		
		uvs[0] = left;
		uvs[1] = right;
		uvs[2] = top;
		
		uvs[3] = new Vector2(left.x,left.y);
		uvs[4] = new Vector2(right.x,right.y);
		uvs[5] = new Vector2(top.x, top.y);
//		
		uvs[6] = left;
		uvs[7] = top;
		uvs[8] = right;
		
		uvs[9] = right;
		uvs[10] = left;
		uvs[11] = top;
		 
		tetraMesh.uv = uvs;
		MeshFilter filter = (MeshFilter)gameObject.GetComponent("MeshFilter");
		filter.mesh = tetraMesh;
		//filter.mesh.RecalculateNormals();
		tetraMesh.RecalculateBounds();
		tetraMesh.Optimize();
	}
	
	void createCube() {
		Mesh cubeMesh = new Mesh();
		Vector3[] cubeVerts = new Vector3[24];
		cubeVerts[0] = new Vector3(0,0,0);
		cubeVerts[1] = new Vector3(1,0,0);
		cubeVerts[2] = new Vector3(0,1,0);
		cubeVerts[3] = new Vector3(1,1,0);
		
		cubeVerts[4] = new Vector3(0,0,1);
		cubeVerts[5] = new Vector3(1,0,1);
		cubeVerts[6] = new Vector3(0,1,1);
		cubeVerts[7] = new Vector3(1,1,1);
		
		cubeVerts[8] = new Vector3(1,0,0);
		cubeVerts[9] = new Vector3(1,1,0);
		cubeVerts[10] = new Vector3(1,1,1);
		cubeVerts[11] = new Vector3(1,0,1);
		
		cubeVerts[12] = new Vector3(0,0,0);
		cubeVerts[13] = new Vector3(0,1,0);
		cubeVerts[14] = new Vector3(0,0,1);
		cubeVerts[15] = new Vector3(0,1,1);
		
		cubeVerts[16] = new Vector3(0,1,0);
		cubeVerts[17] = new Vector3(0,1,1);
		cubeVerts[18] = new Vector3(1,1,0);
		cubeVerts[19] = new Vector3(1,1,1);
		
		cubeVerts[20] = new Vector3(0,0,0);
		cubeVerts[21] = new Vector3(0,0,1);
		cubeVerts[22] = new Vector3(1,0,1);
		cubeVerts[23] = new Vector3(1,0,0);
		cubeMesh.vertices = cubeVerts;
		
		int[] triangles = new int[]{
			0,2,1,2,3,1,
			5,7,4,7,6,4,
			8,9,11,9,10,11,
			14,15,12,15,13,12,
			16,17,18,17,19,18,
			21,20,22,20,23,22
		};
		cubeMesh.triangles = triangles;
		
		Vector2 bottomLeft = new Vector2(0.0f,0.0f);
		Vector2 topLeft = new Vector2(0.0f,1.0f);
		Vector2 topRight = new Vector2(1.0f,1.0f);
		Vector2 bottomRight = new Vector2(1.0f,0.0f);
		Vector2[] uvs = new Vector2[] {
			bottomLeft, bottomRight, topLeft, topRight,
			bottomRight, bottomLeft, topRight, topLeft,
			bottomLeft, topLeft, topRight, bottomRight,
			bottomLeft, bottomRight, topLeft, topRight,
			bottomLeft, bottomRight, topLeft, topRight,
			topLeft, bottomLeft, bottomRight, topRight
		};
		cubeMesh.uv = uvs;
		MeshFilter filter = (MeshFilter)gameObject.GetComponent("MeshFilter");
		filter.mesh = cubeMesh;
		filter.mesh.RecalculateNormals();
	}

//SPHERE
	public void createSphere(float radius, int stacks, int slices) {
		Mesh sphereMesh = new Mesh();
		createSpherePoints(radius, stacks, slices);
		createTriangles(stacks, slices+1);
		createSphereNormals(vertices, radius);
		
		sphereMesh.vertices = vertices;
		sphereMesh.triangles = triangles;
		
		sphereMesh.normals = normals;
		sphereMesh.uv = uvs;
		MeshFilter filter = (MeshFilter)gameObject.GetComponent("MeshFilter");
		filter.mesh = sphereMesh;
	}
	
	private void createSpherePoints(float radius, int stacks, int slices) {
		vertices = new Vector3[(stacks+1)*(slices+1)];
		uvs = new Vector2[vertices.Length];
		float uInc = 1.0f / (slices);
		float vInc = 1.0f / (stacks);
		float thetaInc = 2.0f * Mathf.PI / slices;
		float phiInc = Mathf.PI / stacks;
		for (int i = 0; i <= stacks; i++) {
			float phi = (i * phiInc) - (Mathf.PI / 2.0f);
			float v = i * vInc;
			for (int k = 0; k <= slices; k++) {
				if (k == slices) {
					float theta = 0;
					vertices[i*(slices+1) + k] = sphereToXYZ(theta, phi, radius);
					uvs[i*(slices+1) + k] = new Vector2(1, v);
				} else {
					float theta = k * thetaInc;
					float u = k * uInc;
					vertices[i*(slices+1) + k] = sphereToXYZ(theta, phi, radius);
					uvs[i*(slices+1) + k] = new Vector2(u,v);
				}

			
			}
		}
	}
	
	
	private void createTriangles(int stacks, int slices) {
		triangles = new int[6*stacks*slices];
		for (int st = 0; st < stacks; st++) {
			for (int sl = 0; sl < slices; sl++) {
				int pointIdx = slices*st + sl;
				int nextPointIdx = slices*st + ((sl+1) % slices);
				int upPointIdx = slices*(st+1) + sl;
				int upNextPointIdx = slices*(st+1) + ((sl+1) % slices);
				addToTriangles(new int[]{upPointIdx, nextPointIdx, pointIdx, upPointIdx, upNextPointIdx, nextPointIdx});
			}
		}
	}
	
	private void createSphereNormals(Vector3[] verts, float radius) {
		normals = new Vector3[verts.Length];
		for (int i = 0; i < verts.Length; i++) {
			normals[i] = Vector3.Scale(verts[i], new Vector3(1.0f / radius, 1.0f / radius, 1.0f / radius));
		}
	}
		
	private Vector3 sphereToXYZ(float theta, float phi, float radius) {
		Vector3 point = new Vector3();
		point.x = radius * Mathf.Cos(theta) * Mathf.Cos(phi);
		point.y = radius * Mathf.Sin(phi);
		point.z = radius * Mathf.Sin(theta) * Mathf.Cos(phi);
		return point;
	}
	
	private void addToTriangles(int[] points) {
		for (int i = 0; i < points.Length; i++) {
			triangles[trianglesIndex++] = points[i];
		}
	}
//CYLINDER	
	public void createCylinder(float radius, float height, int slices, GameObject go) {
		Matrix4x4 rot = Matrix4x4.TRS(Vector3.zero, Quaternion.Euler(270f, 0f, 0f), Vector3.one);
		createCylinder(radius, height, slices, go, rot);
	}
	public void createCylinder(float radius, float height, int slices, GameObject go, Matrix4x4 matrix) {
		Mesh cylinderMesh = new Mesh();
		vertices = new Vector3[(slices+1) * 4];
		Vector3[] cylPoints1 = createCylinderPoints(radius, height, slices, true);		
		Vector3[] cylPoints2 = createCylinderPoints(radius, height, slices, false);
		for (int i = 0; i <cylPoints1.Length; i++) {
			vertices[i] = cylPoints1[i];	
		}
		for (int k = 0; k < cylPoints2.Length; k++) {
			vertices[k + cylPoints1.Length] = cylPoints2[k];	
		}
		
		createCylinderNormals(vertices, radius);
		for (int j = 0; j < vertices.Length; j++) {
			vertices[j] = matrix.MultiplyPoint(vertices[j]);
			normals[j] = matrix.MultiplyPoint(normals[j]);
		}
		
		trianglesIndex = 0;
		createCylinderTriangles(slices+1);
		
		cylinderMesh.vertices = vertices;
		cylinderMesh.triangles = triangles;
		cylinderMesh.normals = normals;
		cylinderMesh.uv = uvs;
		
		MeshFilter filter = (MeshFilter)go.GetComponent("MeshFilter");
		filter.mesh = cylinderMesh;
	}
	
	private void createCylinderTriangles(int slices) {
		triangles = new int[12*slices*2];
		for (int p = 0; p < 2; p++) {
			if (p == 0) {
				for (int st = 0; st < 2 - 1; st++) {
					for (int sl = 0; sl < slices; sl++) {
						int pointIdx = slices*st + sl;
						int nextPointIdx = slices*st + ((sl+1) % slices);
						int upPointIdx = slices*(st+1) + sl;
						int upNextPointIdx = slices*(st+1) + ((sl+1) % slices);
						addToTriangles(new int[]{upPointIdx, pointIdx, nextPointIdx, upPointIdx, nextPointIdx, upNextPointIdx});
					}
				}
			} else {
				for (int i = 0; i < 2; i++) {
					bool top = (i == 0) ? false : true;
					for (int k = 0; k < slices - 1; k++) {
						triangles[trianglesIndex++] = i*slices + k + 1;
						if (!top) {
							triangles[trianglesIndex++] = i*slices + k;
							triangles[trianglesIndex++] = i*slices;
						} else {
							triangles[trianglesIndex++] = i*slices;
							triangles[trianglesIndex++] = i*slices + k;
						}
					}
				}	
			}
		}
	}
	
	private Vector3[] createCylinderPoints(float radius, float height, int slices, bool capUV) {
		Vector3[] cylPoints = new Vector3[(slices+1) * 2];
		uvs = new Vector2[vertices.Length];
		float uInc = 1.0f / (slices);
		float vInc = 1;
		float thetaInc = 2.0f * Mathf.PI / slices;
		if (capUV) {
			uInc = 1.0f / (slices / 4.0f);	
			vInc = 1.0f/ (slices / 2.0f);
		}
		for (int i = 0; i < 2; i++) {
			for (int k = 0; k <= slices; k++) {
				if (k == slices) {
					cylPoints[i*(slices+1) + k] = new Vector3(radius * Mathf.Cos(0), i * -height, radius * Mathf.Sin(0));
					uvs[i*(slices+1) + k] = new Vector2(1, i * vInc);
				} else {
					Vector3 point = new Vector3();
					point.x = radius * Mathf.Cos(thetaInc * k);
					point.y = i * -height;
					point.z = radius * Mathf.Sin(thetaInc * k);
					cylPoints[i*(slices+1) + k] = point;
					float u = k * uInc;
					uvs[i*(slices+1) + k] = new Vector2(u, i * vInc);
				}
			}
		}
		return cylPoints;
	}
	
	private void createCylinderNormals(Vector3[] verts, float radius) {
		normals = new Vector3[verts.Length];
		for (int i = 0; i < verts.Length; i++) {
			normals[i] = Vector3.Scale(verts[i], new Vector3(1.0f / radius, 1.0f / radius, 1.0f / radius));
			normals[i].y = 0.0f;
		}
	}
//CAPSULE
	public void createCapsule(float radius, float height, int stacks, int slices) {
		Mesh capsuleMesh = new Mesh();
		createCapsulePoints(radius, stacks, slices, height);
		createCapsuleTriangles(stacks, slices+1);
		createSphereNormals(vertices, radius);
		createCapsuleNormals(vertices, radius, stacks, slices, height);
		
		capsuleMesh.vertices = vertices;
		capsuleMesh.triangles = triangles;
		capsuleMesh.uv = uvs;
		capsuleMesh.normals = normals;
		MeshFilter filter = (MeshFilter)gameObject.GetComponent("MeshFilter");
		filter.mesh = capsuleMesh;

	}
	
	private void createCapsulePoints(float radius, int stacks, int slices, float height) {
		vertices = new Vector3[stacks*(slices+1)];
		uvs = new Vector2[vertices.Length];
		float uInc = 1.0f / (slices);
		float vInc = 1.0f / (stacks);
		float thetaInc = 2.0f * Mathf.PI / slices;
		float phiInc = Mathf.PI / (stacks - 1);
		for (int i = 0; i < stacks; i++) {
			float v = i * vInc;
			float heightJump = 0;
			float vJump = 0;
			if (i > stacks / 2) {
				heightJump = height;
				vJump = 1;
			}
			float phi = i * phiInc - Mathf.PI / 2.0f;
			for (int k = 0; k <= slices; k++) {
				if (k == slices) {
					float theta = 0;
					vertices[i*(slices+1) + k] = sphereToXYZ(theta, phi, radius);
					vertices[i*(slices+1) + k].y += heightJump;
					uvs[i*(slices+1) + k] = new Vector2(1, v + vJump);

				} else {
					float u = k * uInc;
					float theta = k * thetaInc;
					vertices[i*(slices+1) + k] = sphereToXYZ(theta, phi, radius);
					vertices[i*(slices+1) + k].y += heightJump;
					uvs[i*(slices+1) + k] = new Vector2(u,v + vJump);
				}
			}
		}
	}
	
	private void createCapsuleTriangles(int stacks, int slices) {
		triangles = new int[6*stacks*slices];
		for (int st = 0; st < (stacks - 1); st++) {
			for (int sl = 0; sl < slices; sl++) {
				int pointIdx = slices*st + sl;
				int nextPointIdx = slices*st + ((sl+1) % slices);
				int upPointIdx = slices*(st+1) + sl;
				int upNextPointIdx = slices*(st+1) + ((sl+1) % slices);
				addToTriangles(new int[]{upPointIdx, nextPointIdx, pointIdx, upPointIdx, upNextPointIdx, nextPointIdx});
			}
		}
	}
	
	private void createCapsuleNormals(Vector3[] verts, float radius, int stacks, int slices, float heightJump) {
		normals = new Vector3[verts.Length];
		for (int i = 0; i < stacks; i++) {
			bool zeroY = false;
			if ((i == stacks / 2) || (i == (stacks / 2 + 1))) {
				zeroY = true;	
			}
			for (int k = 0; k < slices; k++) {
				Vector3 pointVec = verts[i*slices + k];
				if (i > stacks / 2) {
					pointVec.y -= heightJump;	
				}
				normals[i*slices + k] = Vector3.Scale(pointVec, new Vector3(1.0f / radius, 1.0f / radius, 1.0f / radius));
				if (zeroY == true) {
					normals[i*slices + k].y = 0.0f;
				}
			}
		}
	}
//TAURUS
	public void createTaurus(float smallRadius, float largeRadius, int slices, int points) {
		Mesh taurusMesh = new Mesh();
		vertices = new Vector3[(points+1)*slices];
		uvs = new Vector2[vertices.Length];
		float uInc = 1.0f / (slices);
		float vInc = 1.0f / (points);
		float thetaInc = 2.0f * Mathf.PI / (slices - 1);
		float phiInc = 2.0f * Mathf.PI / points;
		for (int i = 0; i < slices; i++) {
			float theta = i * thetaInc;
			float v = i * vInc;
			for (int k = 0; k <= points; k++) {
				if (k == points) {
					float phi = 0;
					float u = 1;
					vertices[i*(slices+1) + k] = taurusToXYZ(theta, phi, smallRadius, largeRadius);
					uvs[i*(slices+1) + k] = new Vector2(u,v);
				} else {
					float phi = k * phiInc;
					float u = k * uInc;
					vertices[i*(slices+1) + k] = taurusToXYZ(theta, phi, smallRadius, largeRadius);
					uvs[i*(slices+1) + k] = new Vector2(u,v);
				}
			}
		}
		
		createTaurusTriangles(slices,points+1);
		createTaurusNormals(slices, points, largeRadius, smallRadius);
		taurusMesh.vertices = vertices;
		taurusMesh.triangles = triangles;
		taurusMesh.normals = normals;
		taurusMesh.uv = uvs;
		
		MeshFilter filter = (MeshFilter)gameObject.GetComponent("MeshFilter");
		filter.mesh = taurusMesh;;

	}
	
	private void createTaurusNormals(int slices, int points, float largeRadius, float smallRadius) {
		float thetaInc = 2.0f * Mathf.PI / (slices - 1);
		normals = new Vector3[vertices.Length];
		for (int i = 0; i < slices; i++) {
			float theta = i * thetaInc;
			for (int k = 0; k < points; k++) {
				Vector3 distanceFromCenter = new Vector3();
				distanceFromCenter.x = largeRadius * Mathf.Cos(theta);
				distanceFromCenter.y = 0.0f;
				distanceFromCenter.z = largeRadius * Mathf.Sin(theta);
				Vector3 subtractedPoint = vertices[i*slices + k] - distanceFromCenter;
				Vector3 normalizedPoint = Vector3.Scale(subtractedPoint, new Vector3(1.0f / smallRadius, 1.0f / smallRadius, 1.0f / smallRadius));
				normals[i*slices + k] = normalizedPoint;
	
			}
		}
	}
	
	private Vector3 taurusToXYZ(float theta, float phi, float smallRadius, float largeRadius) {
		Vector3 distanceFromCenter = new Vector3();
		distanceFromCenter.x = largeRadius * Mathf.Cos(theta);
		distanceFromCenter.y = 0.0f;
		distanceFromCenter.z = largeRadius * Mathf.Sin(theta);
		Vector3 point = new Vector3();
		point.x = smallRadius * Mathf.Cos(theta) * Mathf.Cos(phi);
		point.y = smallRadius * Mathf.Sin(phi);
		point.z = smallRadius * Mathf.Sin(theta) * Mathf.Cos(phi);
		point += distanceFromCenter;
		return point;	
	}
	
	private void createTaurusTriangles(int stacks, int slices) {
		triangles = new int[6*stacks*slices];
		for (int st = 0; st < stacks - 1; st++) {
			for (int sl = 0; sl < slices; sl++) {
				int pointIdx = slices*st + sl;
				int nextPointIdx = slices*st + ((sl+1) % slices);
				int upPointIdx = slices*(st+1) + sl;
				int upNextPointIdx = slices*(st+1) + ((sl+1) % slices);
				addToTriangles(new int[]{upPointIdx, pointIdx, nextPointIdx, upPointIdx, nextPointIdx, upNextPointIdx});
			}
		}
	}

}

