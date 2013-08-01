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
	public Material mat;

	// Use this for initialization
	void Start () {		
		cylinders = new GameObject[32];
	}
	
	public GameObject drawCylinderWithEndpoints(Vector3 startPoint, Vector3 endPoint) {
		//if (startPoint != endPoint) {
			float xD = endPoint.x - startPoint.x;
			float yD = endPoint.y - startPoint.y;
			float zD = endPoint.z - startPoint.z;
			float distance = Mathf.Sqrt(xD*xD + yD*yD + zD*zD);
			
			/*if (GameObject.Find(System.String.Format("{0}", nextIndex))) {
				cylinderPrefab = GameObject.Find(System.String.Format("{0}", nextIndex));
				cylinderPrefab.transform.position = Vector3.zero;
				cylinderPrefab.transform.rotation = Quaternion.identity;
			} else {*/
				cylinderPrefab = (GameObject)Instantiate(cylinderPrefab, Vector3.zero, Quaternion.identity);
			//}
			createCylinder(0.05f, distance, 8, cylinderPrefab);
			cylinderPrefab.renderer.material = mat;
			cylinderPrefab.transform.Translate(startPoint);
			cylinderPrefab.transform.LookAt(endPoint);	
			cylinderPrefab.name = System.String.Format("{0}", nextIndex);
			//cylinders[nextIndex++] = cylinderPrefab;
			return cylinderPrefab;
		//}			 
	}
		
	float calculateRotation(Vector2 startPoint, Vector2 endPoint) {
		Vector2 legStartPoint = new Vector2(endPoint.x, startPoint.y);
		Vector2 legEndPoint = new Vector2(endPoint.x, endPoint.y);
		float legDistance = Mathf.Sqrt(Mathf.Pow((legEndPoint.x - legStartPoint.x), 2) + Mathf.Pow((legEndPoint.y - legStartPoint.y), 2));
		float hypDistance = Mathf.Sqrt(Mathf.Pow((endPoint.x - startPoint.x), 2) + Mathf.Pow((endPoint.y - startPoint.y), 2));
		
		float angle = Mathf.Asin(legDistance / hypDistance) * (180.0f / Mathf.PI);
		return angle;
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
}

