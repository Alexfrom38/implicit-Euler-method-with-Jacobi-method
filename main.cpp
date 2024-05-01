#include<iostream>
#include <fstream>
#include<cmath>
#include<omp.h>
#include<iomanip>
#include<cmath>

inline void Matrix_Vector_mult(double* elements, size_t* jptr, size_t* iptr, double* vector_result, double* vector_operand, size_t size_of_matrix)
{
	// vectors and arrays already exist
	for (size_t i = 0; i < size_of_matrix; i++)
	{
		vector_result[i] = 0;
		for (size_t j = iptr[i]; j < iptr[i + 1]; j++)
			vector_result[i] = vector_result[i] + elements[j - 1] * vector_operand[jptr[j - 1]];

	}

}

inline void print_iptr(size_t* iptr, size_t size_of_matrix)
{
	for (size_t j = 0; j < size_of_matrix + 1; j++)
		std::cout << iptr[j] << ' ';
	std::cout << std::endl;
}

inline void print_jptr(size_t* jptr, size_t count_elements)
{
	for (size_t j = 0; j < count_elements; j++)
		std::cout << jptr[j] << ' ';
	std::cout << std::endl;
}

inline void print_elements(double* elements, size_t count_elements)
{
	for (size_t j = 0; j < count_elements; j++)
		std::cout << elements[j] << ' ';
	std::cout << std::endl;
}

void Insert_In_File(double* new_array, size_t count, std::fstream& stream)
{
	if (stream.is_open())
	{
		for (size_t i = 0; i < count; i++)
			stream << new_array[i] << std::setprecision(16) << " ";
		stream << "\n";
		stream << "\n";
	}
	else
		throw "The file isn't exist";
}

inline size_t get_count_elements_in_matrix(double* matrix, size_t size_of_matrix)
{
	size_t count_elements = 0;
	for (size_t i = 0; i < size_of_matrix * size_of_matrix; i++)
		if (matrix[i] != 0)
			count_elements++;
	return count_elements;
}

inline void fill_sparce_matrix(double* matrix, size_t size_of_matrix, double* elements, size_t* jptr, size_t* iptr)
{
	int NNZ = 0;
	size_t counter = 0;
	iptr[0] = 1;
	for (size_t i = 0; i < size_of_matrix; i++)
	{
		if (i > 0)
			iptr[i] = NNZ + 1;

		for (size_t j = 0; j < size_of_matrix; j++)
		{
			if (matrix[i * size_of_matrix + j] != 0)
			{
				elements[counter] = matrix[i * size_of_matrix + j];
				jptr[counter] = j;
				counter++;
				NNZ++;
			}
		}
	}
	iptr[size_of_matrix] = NNZ + 1;

}

void Jacobi(double* elements, size_t* iptr, size_t* jptr, size_t size_of_matrix, double* vector_approx, double* column_of_free_numbs, double* TempX, double* norm, double accur)
{
	double norm_max = -1;

	do
	{

		for (size_t i = 0; i < size_of_matrix; i++)
		{
			TempX[i] = column_of_free_numbs[i];
			for (size_t j = iptr[i]; j < iptr[i + 1]; j++)
			{
				if (i != jptr[j - 1])
					TempX[i] -= elements[j - 1] * vector_approx[jptr[j - 1]];
			}
			for (size_t j = iptr[i]; j < iptr[i + 1]; j++)
			{
				if (i == jptr[j - 1])
					TempX[i] = TempX[i] / elements[j - 1];
			}
		}

		for (size_t r = 0; r < size_of_matrix; r++)
			std::cout << TempX[r] << std::setprecision(16) << ' ';
		std::cout << std::endl;

		Matrix_Vector_mult(elements, jptr, iptr, norm, vector_approx, size_of_matrix);

		// norm_max = fabs(vector_approx[0] - TempX[0]);
		norm_max = fabs(norm[0] - column_of_free_numbs[0]);
		for (size_t q = 0; q < size_of_matrix; q++)
		{
			//if (fabs(vector_approx[q] - TempX[q]) > norm_max)
			if (fabs(norm[q] - column_of_free_numbs[q]) > norm_max)
				norm_max = fabs(norm[q] - column_of_free_numbs[q]);

			vector_approx[q] = TempX[q];
		}

	} while (norm_max > accur);

}


int main()
{
	double X_max = 0.0;
	double deltaX = 0.0;
	double T_max = 0.0;
	double deltaT = 0.0;
	double sigma = 0.0;
	double epsilon = 0.000001;
	size_t count = 0;
	double time = 0.0;
	double tmp = 0.0;
	double norm = 0;
	int columns = 0;
	int rows = 0;
	int num_threads = 1;
	size_t count_of_iter = 0;
	double A_i = 0, B_i = 0, C_i = 0;

	const double PI = 3.141592653589793;


	//////////////////////////////////////////////////////////////////////////////////////////
	std::fstream f;
	f.open("output.txt", std::fstream::in | std::fstream::out);
	std::ifstream f_in;
	//////////////////////////////////////////////////////////////////////////////////////////
	f_in.open("const_initial.txt");
	if (f_in.is_open())
	{
		f_in >> X_max;
		f_in >> deltaX;
		f_in >> T_max;
		f_in >> deltaT;
		f_in >> sigma;
		f_in >> epsilon;
		f_in >> num_threads;

	}
	else
		throw "file wasn't opened";
	f_in.close();
	//////////////////////////////////////////////////////////////////////////////////////////   
	A_i = (-1 * sigma * deltaT) / (deltaX * deltaX);
	B_i = 1 + (2 * sigma * deltaT) / (deltaX * deltaX);
	C_i = (-1 * sigma * deltaT) / (deltaX * deltaX);
	//////////////////////////////////////////////////////////////////////////////////////////
	if (X_max == 0.0 || deltaX == 0.0)
		throw "count = 0 or count is infinity";
	//////////////////////////////////////////////////////////////////////////////////////////
	count = static_cast<size_t>((int)(X_max / deltaX)) + 1;
	//////////////////////////////////////////////////////////////////////////////////////////
	double* temp_array = nullptr;
	temp_array = new double[count];

	double* curr_array = nullptr;
	curr_array = new double[count];

	double* tmp_buffer = nullptr;
	tmp_buffer = new double[count];

	double* matrix = nullptr;
	matrix = new double[count * count];
	//////////////////////////////////////////////////////////////////////////////////////////
	double t1 = 0, t2 = 0, dt = 0;
	//////////////////////////////////////////////////////////////////////////////////////////
	for (size_t q = 0; q < count; q++)
		temp_array[q] = sin(PI * deltaX * q);

	temp_array[count - 1] = 0.0;
	Insert_In_File(temp_array, count, f);
	//////////////////////////////////////////////////////////////////////////////////////////

	t1 = omp_get_wtime();
	omp_set_num_threads(num_threads);
#pragma omp parallel for shared (matrix) 
	for (int r = 0; r < count * count; r++)
	{
		matrix[r] = 0;
	}

#pragma omp parallel for shared (matrix) 
	for (int j = 0; j < count; j++)
	{
		if (j - 1 > 0)
			matrix[(j - 1) + j * count] = A_i;
		matrix[j + j * count] = B_i;
		if (j + 1 < count)
			matrix[(j + 1) + j * count] = C_i;

	}

	matrix[0] = 1;
	matrix[count * count - 2] = 0;
	matrix[count * count - 1] = 1;
	matrix[1] = 0;
	matrix[count] = A_i;



#pragma omp parallel for shared (curr_array,temp_array) 
	for (int q = 0; q < count; q++)
		curr_array[q] = temp_array[q];

	while (time <= T_max)
	{

		//////////////////////////////////////////////////////////////////////////////////////////
		do
		{

#pragma omp parallel shared (tmp_buffer,temp_array,matrix,curr_array) 
			{
#pragma	omp for
				for (int i = 0; i < count; i++) {
					tmp_buffer[i] = temp_array[i];
					for (int g = 0; g < count; g++)
					{
						if (i != g)
							tmp_buffer[i] -= matrix[i * count + g] * curr_array[g];
					}
					tmp_buffer[i] /= matrix[i * count + i];
				}
			}
			norm = fabs(curr_array[0] - tmp_buffer[0]);

#pragma omp parallel for shared (curr_array,tmp_buffer) 
			for (int h = 0; h < count; h++)
			{
				if (fabs(curr_array[h] - tmp_buffer[h]) > norm)
					norm = fabs(curr_array[h] - tmp_buffer[h]);

				curr_array[h] = tmp_buffer[h];
			}
			count_of_iter++;
		} while (norm > epsilon);

		std::cout << count_of_iter << std::endl;
		count_of_iter = 0;
		//////////////////////////////////////////////////////////////////////////////////////////
#pragma omp parallel for shared (temp_array,curr_array) 
		for (int mm = 0; mm < count; mm++)
			temp_array[mm] = curr_array[mm];


		time = time + deltaT;
	}



	t2 = omp_get_wtime();
	dt = t2 - t1;
	std::cout << dt << std::endl;

	Insert_In_File(curr_array, count, f);
	f.close();
	delete[] temp_array;
	delete[] matrix;
	delete[] curr_array;


	return 0;
}

