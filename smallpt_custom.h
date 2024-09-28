#pragma once
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "erand48.h"

// 原版smallpt有一些变量使用了但没在原版代码里定义
#define M_PI 3.141592653589793238462643

// 自定义Smallpt主类
class SmallptCustom
{
public:
	void Execute();
};

// Vector结构体 用于颜色、三维点、向量等数据量
struct Vector
{
	double x, y, z;

	Vector(double xValue = 0, double yValue = 0, double zValue = 0)
	{
		x = xValue;
		y = yValue;
		z = zValue;
	}

	Vector operator+(const Vector& outVec) const { return Vector(x + outVec.x, y + outVec.y, z + outVec.z); }
	Vector operator-(const Vector& outVec) const { return Vector(x - outVec.x, y - outVec.y, z - outVec.z); }
	Vector operator*(double k) const { return Vector(x * k, y * k, z * k); }
	Vector mult(const Vector& outVec) const { return Vector(x * outVec.x, y * outVec.y, z * outVec.z); }
	Vector& norm() { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }
	double dot(const Vector& outVec) const { return x * outVec.x + y * outVec.y + z * outVec.z; }
	// 原版smalltpt叉乘用%重载 为了简化代码量 但我这里就用cross
	Vector cross(Vector& outVec) { return Vector(y * outVec.z - z * outVec.y, z * outVec.x - x * outVec.z, x * outVec.y - y * outVec.x); }
};

// 三种材质定义 其实我觉得这里不完全算材质 应该算会发生高光反射 折射反射同时发生 漫反射情形的定义
enum Refl_t
{
	DIFF,
	SPEC,
	REFR
};

// 光线定义
struct Ray
{
	Vector origin, dir;
	Ray(Vector originValue, Vector dirValue) :
		origin(originValue),
		dir(dirValue)
	{}
};

// 场景中几何体很简单 一堆球体 球体和光线求交推导也比较简单 下面是球体数据结构定义
struct Sphere
{
	double radius; // 半径
	Vector pos; // 球心位置

	Vector emission, color; // 表面自发光强度 表面散射颜色(物体固有颜色 光泽反射时颜色不体现 但漫反射会改变入射光的颜色)

	Refl_t refl;  // 对光的反射类型

	Sphere(double radiusValue, Vector posValue, Vector emissionValue, Vector colorValue, Refl_t reflValue):
		radius(radiusValue),
		pos(posValue),
		emission(emissionValue),
		color(colorValue),
		refl(reflValue)
	{
	}

	// 返回传入的光线与该球体接触的距离（光线是一个向量 距离就够算出真实接触位置了）返回0则表示未接触
	double intersect(const Ray& ray) const 
	{
		// 实际就是将光线向量代入球体方程 通过求根公式求解
		// (Dir*Dir)t^2 + 2Dir dot (origin - pos)t + (origin - pos) dot (origin - pos) - r^2 = 0
		Vector op = pos - ray.origin;  // 这里实际算反了 是po向量 但最后返回值求根公式也取反了 所以结论是一样的
		double eps = 1e-4; // 数学中表明足够小的数
		double intersectDis = 1e-4;
		double b = op.dot(ray.dir);  // 求根公式中的b
		double det = b * b - op.dot(op) + radius * radius;
		if (det < 0)
		{
			return 0;
		}
		else
		{
			det = sqrt(det);
			return (intersectDis = b - det) > eps ? intersectDis : ((intersectDis = b + det) > eps ? intersectDis : 0);  // return smallpt positive intersectDis
		}
	}
};

// 光线穿过的虚拟成像屏幕
struct Image
{
	int w, h;  // image size
	Ray cam;  // 相机

	// 基本上通过这几个变量就可以定位相机的方向 以及光线射出后在成像屏幕击中的点对应的局部坐标系用 世界坐标的表示
	Vector cx;  // 就是相机水平方向
	Vector cy;  // 相机的竖直方向
	Vector* targetColor;
	int numSpheres = 9;  // 球体数量
	
	Sphere spheres[9] = {
		// 构成墙壁的六个球体 足够大的球体在局部区域看似平面
		Sphere(1e5,	Vector(1e5 + 1,40.8,81.6),		Vector(), Vector(.75,.25,.25),	DIFF),//Left
		Sphere(1e5,	Vector(-1e5 + 99,40.8,81.6),	Vector(), Vector(.25,.25,.75),	DIFF),//Right
		Sphere(1e5, Vector(50,40.8,1e5),			Vector(), Vector(.75,.75,.75),	DIFF),//Back
		Sphere(1e5, Vector(50,40.8,-1e5 + 170),		Vector(), Vector(),				DIFF),//Front
		Sphere(1e5, Vector(50,1e5,81.6),			Vector(), Vector(.75,.75,.75),	DIFF),//Bottom
		Sphere(1e5, Vector(50,-1e5 + 81.6,81.6),	Vector(), Vector(.75,.75,.75),	DIFF),//Top
		// 场景里的球体
		Sphere(16.5, Vector(27,16.5,47),	Vector(),	Vector(1,1,1) * .999,		SPEC),//Mirr
		Sphere(16.5, Vector(73,16.5,78),	Vector(),	Vector(1,1,1) * .999,		REFR),//Glas
		// 光源
		Sphere(600,	Vector(50,681.6 - 0.27,81.6),	Vector(12,12,12), Vector(),    DIFF),//Lite
	};

	// 单个像素超采样量级
	int samps = 1;
	Image(int w_, int h_, Ray cam_, int samps_):
		w(w_),
		h(h_),
		cam(cam_),
		samps(samps_)
	{
		// 0.5135决定相机视野 smallpt作者估计基于预期效果大概调试出的值
		cx = Vector(w * .5135 / h);
		cy = cx.cross(cam.dir).norm() * 0.5135;
		targetColor = new Vector[w * h];
	}

	inline double clamp(double x)
	{
		return x < 0 ? 0 : (x > 1 ? 1 : x);
	}

	inline int toInt(double x)
	{
		// gamma校正 现阶段成像设备应该不需要均衡了 那目的就是为了显示更符合人眼规律？
		return int(pow(clamp(x), 1 / 2.2) * 255 + .5);
	}

	void outputImage()
	{
#pragma omp parallel for schedule(dynamic, 1) private(r) // OpenMP

		Vector r;  // 存储单个子像素采样结果

		for (int y = 0; y < h; y++)  // 按列遍历所有像素
		{
			fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100. * y / (h - 1));
			unsigned short Xi[3] = { 0, 0, y * y * y };
			for (int x = 0; x < w; x++)  // 遍历某一行像素
			{
				for (int sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++)
				{
					for (int sx = 0; sx < 2; sx++, r = Vector()) 
					{
						for (int s = 0; s < samps; s++)
						{
							double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
							double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
							Vector lightDir = cx * (((sx + 0.5 + dx) / 2 + x) / w - 0.5) + cy * (((sy + 0.5 + dy) / 2 + y) / h - 0.5) + cam.dir;
							Ray outRay = Ray(cam.origin + lightDir * 140, lightDir.norm());
							r = r + radiance(outRay, 0, Xi) * (1. / samps);  // 取单个子像素多次采样后的平均值
						}
						targetColor[i] = targetColor[i] + Vector(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
					}
				}
			}
		}

		// 结果写入ppm后缀文件 方便查看结果
		FILE* f;
		errno_t err = fopen_s(&f, "image.ppm", "w"); // Write image to PPM file.
		fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
		for (int i = 0; i < w * h; i++)
		{
			fprintf(f, "%d %d %d ", toInt(targetColor[i].x), toInt(targetColor[i].y), toInt(targetColor[i].z));
		}
	}

	inline bool intersect(const Ray& ray, double& t, int& id)
	{
		double d, inf = t = 1e20;
		double n = sizeof(spheres) / sizeof(Sphere);

		// 遍历所有球体
		for (int i = int(n); i--;)
		{
			if ((d = spheres[i].intersect(ray)) && d < t)
			{
				t = d;
				id = i;
			}
		}

		return t < inf;
	}

	Vector radiance(const Ray& ray, int depth, unsigned short* Xi, int E = 1)
	{

		double t;   // distance to intersection
		int id = 0; // id of intersected object

		if (!intersect(ray, t, id))
		{
			return Vector(); // 未相交
		}

		const Sphere& obj = spheres[id]; // the hit object
		if (depth > 10)
		{
			return obj.emission;
		}

		Vector x = ray.origin + ray.dir * t;  // 交点位置
		Vector n = (x - obj.pos).norm();  // 法线
		Vector nl = n.dot(ray.dir) < 0 ? n : n * -1;  // 主要用于判断射线是否是在球体内部
		Vector objColor = obj.color;

		double p = objColor.x > objColor.y && objColor.x > objColor.z ? objColor.x : objColor.y > objColor.z ? objColor.y : objColor.z; // max refl

		if (++depth > 5)
		{
			if (erand48(Xi) < p)
				objColor = objColor * (1 / p);
			else
				return obj.emission;
		}

		if (obj.refl == DIFF)
		{  // 理想漫反射

			// 余弦重要性采样
			double r1 = 2 * M_PI * erand48(Xi), r2 = erand48(Xi), r2s = sqrt(r2);
			Vector w = nl, u = ((fabs(w.x) > .1 ? Vector(0, 1) : Vector(1)).cross(w)).norm(), v = w.cross(u);
			Vector d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
			return obj.emission + objColor.mult(radiance(Ray(x, d), depth, Xi));
		}
		else if (obj.refl == SPEC)
		{  // 理想镜面反射
			Ray newOutRay(x, (ray.dir - n * 2 * n.dot(ray.dir)));
			return obj.emission + objColor.mult(radiance(newOutRay, depth, Xi));
		}
		else
		{
			Ray reflRay(x, (ray.dir - n * 2 * n.dot(ray.dir)));
			bool into = n.dot(nl) > 0;
			double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = ray.dir.dot(nl), cos2t;
			if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)  // 只在球体内全反射
				return obj.emission + objColor.mult(radiance(reflRay, depth, Xi));

			Vector tdir = (ray.dir * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();  //折射光线归一化矢量 tdir
			/*
			* 对玻璃材质，前几次反射、折射对最终的结果的贡献比较大，所以选择都进行追踪。但是这样做一下子使得某些追踪的采样次数发散为了原来的
			* 4 倍，如果深度改为 4，可能就是 16 倍。所以代码把发散深度限制在 2 这样一个比较小的数值，也算是在效果和性能间达到一个平衡
			*/
			double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : tdir.dot(n));
			double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
			return obj.emission + objColor.mult(depth > 2 ? (erand48(Xi) < P ?   // Russian roulette
				radiance(reflRay, depth, Xi) * RP : radiance(Ray(x, tdir), depth, Xi) * TP) :
				radiance(reflRay, depth, Xi) * Re + radiance(Ray(x, tdir), depth, Xi) * Tr);
		}
	}
};

