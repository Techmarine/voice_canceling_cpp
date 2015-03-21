#include <string>
#include <fstream>
#include <vector>
using namespace std;

const int n = 14; //number of fractional bits

//This function writes matrix to .csv file with specified filename
void Matrix2CSV(string filename, vector<vector<double>> matrix)
{
	unsigned long rows = matrix.size(); //number of rows
	unsigned long cols = matrix[0].size(); //number of columns

	ofstream file;
	file.open(filename);
	//Values in a row are separated with commas. The last value is followed by 'new line' symbol
	for (unsigned long i = 0; i < rows; i++)
	{
		for (unsigned long j = 0; j < cols; j++)
		{
			if (j < cols - 1)
				file << matrix[i][j] << ',';
			else
				file << matrix[i][j] << '\n';
		}
	}
	file.close();
}

//This function reads .csv file with specified name to matrix
vector<vector<double>> CSV2Matrix(string filename)
{
	ifstream file;
	string line;
	unsigned long cols = 1; //number of columns (initialization)
	unsigned long rows = 0; //number of rows (initialization)

	file.open(filename);
	//Firstly, count the number of rows and columns
	while (getline(file, line)) //read line by line
	{
		if (rows == 0)
		{
			for (int i = 0; i < line.length(); i++)
			{
				if (line[i] == ',') //count commas
				{
					cols = cols + 1;
				}
			}
		}
		rows = rows + 1; //and count the number of lines passed
	}

	//Reset EOF flag and return to the beginning of file
	file.clear();
	file.seekg(0, ios::beg);

	//Create and initialize matrix variable
	vector<vector<double>> matrix(rows);
	for (unsigned long i = 0; i < rows; i++)
	{
		matrix[i].resize(cols);
	}

	//Parse values to the matrix
	unsigned long m = 0; //column index
	unsigned long n = 0; //row index
	int previous = -1; //previous location of comma (for substr)
	while (getline(file, line)) //read line by line
	{
		for (int k = 0; k < line.length(); k++) //interate through the line
		{
			if (line[k] == ',') //find comma
			{
				matrix[n][m] = stod(line.substr(previous + 1, k - previous)); //convert to double the string between two commas
				m++;
				previous = k; //and memorize the comma's position
			}
		}
		matrix[n][m] = stod(line.substr(previous + 1, line.length() - previous)); //the last value in a row is not processed by the loop above, so it is done here

		n++;
		m = 0;
		previous = -1;
	}
	file.close();

	return matrix;
}

//FIR direct form structure filter implementation
vector<double> FIRDirectFormFilter(vector<double> input, vector<double> numerator)
{
	vector<double> output;
	output.resize(input.size());

	vector<double> buffer; //this buffer contatins all delayed samples of the input signal
	buffer.resize(numerator.size());

	double sum = 0; //auxiliary variable for summation

	for (int i = 0; i < input.size(); i++)
	{		
		//Reset sum variable
		sum = 0;

		//Shift values in buffer
		for (int j = numerator.size() - 1; j > 0; j--)
		{
			buffer[j] = buffer[j - 1];
		}
		buffer[0] = input[i];

		//Multiply and sum
		for (int j = 0; j < numerator.size(); j++)
		{
			sum = sum + buffer[j] * numerator[j];
		}
		output[i] = sum;
	}

	return output;
}

//This function returns the whole column of a matrix
vector<double> GetColumn(vector<vector<double>> matrix, int column_index)
{
	vector<double> vec;
	vec.resize(matrix.size());
	for (int i = 0; i < matrix.size(); i++)
	{
		vec[i] = matrix[i][column_index];
	}
	return vec;
}

//This function combines two vectors into a single matrix
template <class type>
vector<vector<type>> CombineChannels(vector<type> left, vector<type> right)
{
	vector<vector<type>> matrix(left.size());

	for (int i = 0; i < left.size(); i++)
	{
		matrix[i].resize(2);
		matrix[i][0] = left[i];
		matrix[i][1] = right[i];
	}
	return matrix;
}

//This function converts floating point vector to fixed point (modeled as integer) vector
vector<int> Float2FixedVector(vector<double> vec, int n)
{
	vector<int> output(vec.begin(),vec.end());
	output.resize(vec.size());

	for (int i = 0; i < vec.size(); i++)
	{
		output[i] = round(vec[i] * pow(2, n));
		
	}
	return output;
}
//FIR direct form structure filter implementation (fixed point version)
vector<int> FIRDirectFormFilterFixed(vector<int> input, vector<int> numerator, int n)
{
	vector<int> output;
	output.resize(input.size());

	vector<int> buffer; //this buffer contatins all delayed samples of the input signal
	buffer.resize(numerator.size());

	int sum = 0; //auxiliary variable for summation

	for (int i = 0; i < input.size(); i++)
	{
		//Reset sum variable
		sum = 0;

		//Shift values in buffer
		for (int j = numerator.size() - 1; j > 0; j--)
		{
			buffer[j] = buffer[j - 1];
		}
		buffer[0] = input[i];

		//Multiply and sum
		for (int j = 0; j < numerator.size(); j++)
		{
			sum = sum + buffer[j] * numerator[j] / pow(2, n); 
			//Dividing by 2^n here is equivalent to truncation of values in register. 
			//The size of a register that stores a product should be 2*n, but setting it to n will result in truncation of n bits
		}
		output[i] = sum;
	}

	return output;
}

//This function converts fixed point (modeled as integer) matrix to floating point matrix
vector<vector<double>> Fixed2FloatMatrix(vector<vector<int>> matrix,int n)
{
	int rows = matrix.size(); //number of rows
	int cols = matrix[0].size(); //number of columns

	vector<vector<double>> output(rows);
	for (int i = 0; i < rows; i++)
	{
		output[i].resize(cols);
	}
	
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			output[i][j] = matrix[i][j] / pow(2, n);
		}
	}
	return output;
}

int main()
{
	//Read wavedata from file
	vector<vector<double>> wavedata = CSV2Matrix("wavedata.csv");

	//Read filter coefficients from file
	vector<vector<double>> numerator = CSV2Matrix("numerator.csv");

	//Floating point representation//
	//Filter both channels
	vector<double> fdata_left = FIRDirectFormFilter(GetColumn(wavedata, 0), numerator[0]);
	vector<double> fdata_right = FIRDirectFormFilter(GetColumn(wavedata, 1), numerator[0]);

	//Combine channels into single matrix and write it to file
	Matrix2CSV("filtered_data_float.csv", CombineChannels(fdata_left, fdata_right));

	//Fixed point representation
	//Convert floating point to fixed point
	vector<int> wavedata_left_fixed = Float2FixedVector(GetColumn(wavedata, 0), n);
	vector<int> wavedata_right_fixed = Float2FixedVector(GetColumn(wavedata, 1), n);
	vector<int> numerator_fixed = Float2FixedVector(numerator[0], n);

	//Filter both channels
	vector<int> fdata_left_fixed = FIRDirectFormFilterFixed(wavedata_left_fixed, numerator_fixed,n);
	vector<int> fdata_right_fixed = FIRDirectFormFilterFixed(wavedata_right_fixed, numerator_fixed,n);

	//Combine channels into single matrix, convert it back to floating point and write it to file
	Matrix2CSV("filtered_data_fixed.csv", Fixed2FloatMatrix(CombineChannels(fdata_left_fixed, fdata_right_fixed), n));
}

