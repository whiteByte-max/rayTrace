#include "smallpt_custom.h"
#include <iostream>

void SmallptCustom::Execute()
{
	std::cout << "test main" << std::endl;
	Ray cam(Vector(50, 52, 295.6), Vector(0, -0.042612, -1).norm());
	int samps = 50;  // ���������ز�������
	Image* image = new Image(1024, 768, cam, samps);

	// ִ��������Ⱦͼ��
	image->outputImage();
}