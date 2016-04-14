
template <class Data> class Matrix3D
{
public:
	Matrix3D();
	Matrix3D(int dim1, int dim2, int dim3);


	~Matrix3D();

private:
	Data *linear;
};