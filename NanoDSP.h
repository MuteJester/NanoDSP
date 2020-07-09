#pragma once
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>
#include <utility>
#include <random>
#include <limits>
#include <exception>
#include <assert.h>
#include <functional>
#include <fstream>
#define DSP_RANDOM_SIGNAL 9912
#define DSP_EMPTY 9913


typedef std::pair<double, double> DoublePair;
typedef std::pair<int, int> IntPair;
typedef std::pair<int, double> IntDoublePair;


//misc 
double Vector_Median(std::vector<double> num) {
	std::sort(num.begin(), num.end());
	if (num.size() % 2 == 0) {
		return (num[(int)std::round(num.size() / 2)] + num[(int)std::round(num.size() / 2) - 1]) / 2;
	}
	else {

		return num[(int)std::round(num.size() / 2)];
	}
}
double Vector_Mean(std::vector<double> const &num) {
	double sum = 0;
	for (auto i : num) {
		sum += i;
	}
	return sum / num.size();
}

#ifndef RANDOMUTILITISDEFINE
#define RANDOMUTILITISDEFINE

class Random_Utilitis {
public:
	static int Random_INT(int minimum_value, int maximum_value) {
		static std::random_device seed;
		static std::mt19937 random_number_generator(seed());
		std::uniform_int_distribution<size_t> indices(minimum_value, maximum_value - 1);
		return (int)indices(random_number_generator);
	}

	static int Random_INT(int minimum_value, int maximum_value, int _seed) {
		static std::random_device seed;
		static std::mt19937 random_number_generator(_seed);
		std::uniform_int_distribution<size_t> indices(minimum_value, maximum_value - 1);
		return (int)indices(random_number_generator);
	}
	static double Random_DOUBLE(double minimum_value, double maximum_value) {
		static std::random_device seed;
		static std::mt19937 random_number_generator(seed());
		std::uniform_real_distribution<> indices(minimum_value, maximum_value);
		return indices(random_number_generator);
	}
	static double Random_DOUBLE(double minimum_value, double maximum_value, double _seed) {
		static std::random_device seed;
		static std::mt19937 random_number_generator(_seed);
		std::uniform_real_distribution<> indices(minimum_value, maximum_value);
		return indices(random_number_generator);
	}

};


#endif // !RANDOMUTILITISDEFINE
#ifndef  SIPLCOMPLEXSYSTEMS
#define SIPLCOMPLEXSYSTEMS


class Complex {
public:
	double Real;
	double Imaginary;
	Complex() {
		this->Real = 0;
		this->Imaginary = 0;
	}
	Complex(double const &real, double const &imaginary) {
		this->Real = real;
		this->Imaginary = imaginary;
	}
	Complex(const Complex &Copy) {
		this->Real = Copy.Real;
		this->Imaginary = Copy.Imaginary;
	}
	static Complex I() {
		return Complex(0, 1);
	}
	static Complex Zero() {
		return Complex(0, 0);
	}

	void operator=(Complex const &source) {
		this->Real = source.Real;
		this->Imaginary = source.Imaginary;
	}
	void operator=(Complex &source) {
		this->Real = source.Real;
		this->Imaginary = source.Imaginary;
	}

	Complex operator+(Complex const &b) {
		return Complex(this->Real + b.Real, this->Imaginary + b.Imaginary);
	}
	Complex operator-(Complex const &b) {
		return Complex(this->Real - b.Real, this->Imaginary - b.Imaginary);

	}
	Complex operator*(Complex const &b) {
		return Complex(this->Real * b.Real - this->Imaginary * b.Imaginary,
			this->Real * b.Imaginary + this->Imaginary * b.Real);
	}
	Complex operator*(double const &b) {
		return Complex(this->Real * b, this->Imaginary*b);
	}
	Complex operator^(int const &B) {
		for (int i = 1; i < B; i++) {
			*this * *this;
		}
	}
	Complex operator^(Complex x) {
		return this->Log()*(x).Exp();

	}
	bool operator==(Complex const &b) {
		if (b.Imaginary == this->Imaginary && this->Real == b.Real) {
			return true;
		}
		else {
			return false;
		}
	}
	bool operator!=(Complex const &b) {
		return *this == b ? false : true;
	}

	Complex Exp() {
		double expReal = std::exp(Real);
		return  Complex(expReal *  std::cos(Imaginary), expReal * std::sin(Imaginary));
	}
	Complex operator/(Complex const &divisor) {

		if ((divisor.Real == 0.0 && divisor.Imaginary == 0.0)) {
			std::abort();
		}


		if (std::abs(divisor.Real) < std::abs(divisor.Imaginary)) {
			double q = divisor.Real / divisor.Imaginary;
			double denominator = divisor.Real * q + divisor.Imaginary;
			return Complex((this->Real * q + this->Imaginary) / denominator, (this->Imaginary * q - this->Real) / denominator);
		}
		else {
			double q = divisor.Imaginary / divisor.Real;
			double denominator = divisor.Imaginary * q + divisor.Real;
			return Complex((this->Imaginary * q + this->Real) / denominator, (this->Imaginary - this->Real * q) / denominator);

		}
	}
	void Conjugate() {
		this->Imaginary *= -1;
	}
	double Abs() {
		if (std::abs(Real) < std::abs(Imaginary)) {
			if (Imaginary == 0.0) {
				return std::abs(Real);
			}
			double q = Real / Imaginary;
			return std::abs(Imaginary) * std::sqrt(1 + q * q);
		}
		else {
			if (Real == 0.0) {
				return std::abs(Imaginary);
			}
			double q = Imaginary / Real;
			return std::abs(Real) * std::sqrt(1 + q * q);
		}
	}
	Complex Log() {
		return Complex(std::log(Abs()),
			std::atan2(Imaginary, Real));
	}
	Complex Sqrt() {
		if (this->Real == 0.0 && Imaginary == 0.0) {
			return Complex(0.0, 0.0);
		}

		double t = std::sqrt((std::abs(Real) + this->Abs()) / 2.0);
		if (Real >= 0.0) {
			return Complex(t, Imaginary / (2.0 * t));
		}
		else {
			return Complex(std::abs(Imaginary) / (2.0 * t),
				std::copysign(1.0
					, Imaginary) * t);
		}
	}


};


std::vector<Complex> Discrete_Fourier_Transform(std::vector<int> &dat) {
	std::vector<Complex> Res;
	Res.reserve(dat.size());
	int d_size = dat.size();


	for (int k = 0; k < dat.size(); k++) {
		Complex sum(0, 0);
		for (int n = 0; n < dat.size(); n++) {
			Complex exp_power(0, ((-2 * std::_Pi) / d_size) *(k*n));
			exp_power = exp_power.Exp();
			exp_power.Real *= dat[n];
			sum = sum + (exp_power);
		}
		Res.emplace_back(sum);
	}
	return Res;
}
std::vector<Complex> Discrete_Fourier_Transform(std::vector<Complex> &dat) {
	std::vector<Complex> Res;
	Res.reserve(dat.size());
	int d_size = dat.size();
	for (int k = 0; k < dat.size(); k++) {
		Complex sum(0, 0);
		//sum = std::accumulate(dat.begin(), dat.end(), Complex::Zero(), [](Complex &v) {return Complex(0, 0); });

		for (int n = 0; n < dat.size(); n++) {
			Complex exp_power(0, ((-2 * std::_Pi) / d_size) *(k*n));
			exp_power = exp_power.Exp();

			sum = sum + (dat[n] * exp_power);
		}
		Res.emplace_back(sum);
	}
	return Res;
}

std::vector<double> Reverse_Discrete_Fourier_Transform(std::vector<Complex> &dat) {
	std::vector<double> Res;
	Res.reserve(dat.size());
	int d_size = dat.size();
	for (int k = 0; k < dat.size(); k++) {
		Complex sum(0, 0);
		for (int n = 0; n < dat.size(); n++) {
			Complex exp_power(0, ((2 * std::_Pi)*(k*n)) / d_size);
			exp_power = exp_power.Exp();
			sum = sum + (dat[n] * exp_power);
		}
		sum = sum * (Complex(1, 0) / Complex(dat.size(), 0));
		Res.emplace_back(sum.Real);
	}
	return Res;
}
std::vector<Complex> Reverse_Discrete_Fourier_Transform_Complex(std::vector<Complex> &dat) {
	std::vector<Complex> Res;
	Res.reserve(dat.size());
	int d_size = dat.size();
	for (int k = 0; k < dat.size(); k++) {
		Complex sum(0, 0);
		for (int n = 0; n < dat.size(); n++) {
			Complex exp_power(0, ((2 * std::_Pi)*k*n) / d_size);
			exp_power = exp_power.Exp();
			sum = sum + (dat[n] * exp_power);
		}
		sum = sum * (Complex(1, 0) / Complex(dat.size(), 0));
		Res.push_back(sum);
	}
	return Res;
}



int bitReverse(int n, int bits) {
	int reversedN = n;
	int count = bits - 1;

	n >>= 1;
	while (n > 0) {
		reversedN = (reversedN << 1) | (n & 1);
		count--;
		n >>= 1;
	}

	return ((reversedN << count) & ((1 << bits) - 1));
}
void FFT(std::vector<Complex> &Values) {
	int length = Values.size();
	int bits = (int)(std::log(length) / std::log(2));
	for (int j = 1; j < length / 2; j++) {

		int swapPos = bitReverse(j, bits);
		Complex temp = Values[j];
		Values[j] = Values[swapPos];
		Values[swapPos] = temp;
	}

	for (int N = 2; N <= length; N <<= 1) {
		for (int i = 0; i < length; i += N) {
			for (int k = 0; k < N / 2; k++) {

				int evenIndex = i + k;
				int oddIndex = i + k + (N / 2);
				Complex even = Values[evenIndex];
				Complex odd = Values[oddIndex];

				double term = (-2 * std::_Pi * k) / (double)N;

				Complex exp = Complex(std::cos(term), std::sin(term))*(odd);

				Values[evenIndex] = even + (exp);
				Values[oddIndex] = even - (exp);
			}
		}
	}
}

#endif // ! SIPLCOMPLEXSYSTEMS




class Signal {
public:
	std::vector<IntDoublePair> Sig;
	int number_of_samples;
	int time_period=1;
	int start_time=0;

	Signal(){}
	Signal(int number_of_samples,int start_time ,int time_period, int mode = DSP_EMPTY) {
		this->start_time = start_time;
		if (mode == DSP_EMPTY) {
			this->number_of_samples = number_of_samples;
			this->time_period = time_period;
			for (int i = 0; i < number_of_samples; i++) {
				Sig.push_back(IntDoublePair(start_time + i*time_period, 0));
			}
		}
		else if (mode == DSP_RANDOM_SIGNAL) {
			this->number_of_samples = number_of_samples;
			this->time_period = time_period;
			Random_Utilitis rnd;
			for (int i = 0; i < number_of_samples; i++) {
				Sig.push_back(std::pair<int,double>(start_time + i * time_period,rnd.Random_DOUBLE(0,40)));
				//std::cout << start_time + i * time_period << std::endl;
			}
		}
	
	}
	Signal(int number_of_samples, int start_time,int time_period, int mode,int random_seed) {
			Random_Utilitis rnd;
			this->start_time = start_time;

			this->number_of_samples = number_of_samples;
			this->time_period = time_period;
			for (int i = 0; i < number_of_samples; i++) {
				Sig.push_back(IntDoublePair(start_time + i*time_period, rnd.Random_DOUBLE(0, 30, random_seed)));
			}
		
	}
	Signal(std::vector<IntDoublePair> const &signal) {
		this->Sig = signal;
		number_of_samples = (int)signal.size();
		this->time_period = (int)signal[1].first;
	}
	Signal(Signal const &signal) {
		this->Sig = signal.Sig;
		number_of_samples = signal.number_of_samples;
		this->time_period = signal.time_period;
		this->start_time = signal.start_time;
	}
	Signal(std::vector<double> const &data) {
		this->number_of_samples = (int)data.size();
		for (int i = 0; i < data.size(); i++) {
			this->Sig.push_back(IntDoublePair(i, data[i]));
		}
	}
	static Signal Sinusoid(double const &Amplitude, double const &Frequency, double const &Phase,int const &number_of_samples,int const &sample_rate,int const &start_time) {
		Signal Res;
		for (int i = start_time; i < start_time + sample_rate * number_of_samples; i += sample_rate) {
			Res.Sig.push_back(IntDoublePair(i, Amplitude*std::sin(Frequency*(i / sample_rate) + Phase)));
		}
		return Res;
	}
	double &operator[](int at_time_x) {
		for (int i = 0; i < number_of_samples; i++) {
			if (this->Sig[i].first == at_time_x) {
				return this->Sig[i].second;
			}
		}
		std::exception e("Invalid Time");
	}
	void operator()(double const& new_point,int const &at_time) {
		this->number_of_samples++;
		assert(at_time % time_period == 0);
		Sig.push_back(IntDoublePair(at_time, new_point));
	}

	Signal Flip() {
		Signal fliped(*this);
		fliped.start_time = -(fliped.start_time + (fliped.number_of_samples-1)*fliped.time_period);
		for (int i = 0; i < number_of_samples; i++) {
			//fliped.Sig[i].second *= -1;
			fliped.Sig[i].first *= -1;

		}
		return fliped;
	}
	Signal Scale(double const &scale) {
		if (scale < 1) {
			int multiplyer = (1 / scale);

			Signal scaled(number_of_samples*multiplyer, this->start_time*multiplyer, this->time_period);

			for (int i = start_time; i < start_time + number_of_samples * time_period; i += time_period) {
				if ((int)(i / scale) % time_period == 0) {
					scaled[i / scale] = this->operator[](i);
				}
			}
			return scaled;

			
			
		}
		else if (scale >= 1) {
			double ns;
			for (int i = start_time; i < start_time + number_of_samples * time_period; i+=time_period) {
				if (i % (int)scale == 0) {
					ns = i / scale;
					break;
				}
			}
		
			Signal scaled(std::ceil(number_of_samples/scale),ns,this->time_period);

			for (int i = start_time; i < start_time + number_of_samples*time_period ; i+= time_period) {
				if ((int)(i / scale)%time_period ==0 ) {
					scaled[i / scale] = this->operator[](i);
				}
			}
			return scaled;
		}
	}
	Signal Shift(int const &units) {
		Signal Shifted(*this);
		Shifted.start_time += units * time_period;
		for (int i = 0; i < this->number_of_samples; i++) {
			Shifted.Sig[i].first += units*this->time_period;
		}
		return Shifted;
	}
	Signal Get_Even_Part() {
		Signal Even(this->number_of_samples, this->start_time, this->time_period);
		Signal flip = this->Flip();
		for (int i = start_time, j=flip.start_time; i < start_time + time_period * number_of_samples; i += time_period,j+= time_period) {
				Even[i] = 0.5 * (this->operator[](i) + flip[j]);
		}
		return Even;
	}
	Signal Get_Odd_Part() {
		Signal Even(this->number_of_samples, this->start_time, this->time_period);
		Signal flip(*this);
		flip = flip.Flip();
		for (int i = start_time, j = flip.start_time; i < start_time + time_period * number_of_samples; i += time_period, j += time_period) {
			Even[i] = 0.5 * (this->operator[](i) - flip[j]);
		}
		return Even;
	}
	
};



Signal operator+(Signal sig_a, Signal sig_b) {
	assert(sig_b.time_period == sig_a.time_period);
	Signal result(sig_a);

	for (int i = sig_b.start_time; i < sig_b.start_time + sig_b.number_of_samples*sig_b.time_period; i+=sig_b.time_period) {
		if (i < sig_a.start_time || i > (sig_a.start_time +sig_a.number_of_samples*sig_a.time_period)) {
			result(sig_b[i], i);
		}
		else {
			result[i] += sig_b[i];
		}
	}

	return result;

}
std::ostream &operator<<(std::ostream &out, IntDoublePair const &p) {
	out << "[Time:  " << p.first << ", Value: " << p.second << " ] ";
	return out;
}


std::ostream &operator<<(std::ostream &out,Signal  &sig) {
	for (int i = 0; i < sig.number_of_samples; i++) {
		//out << sig[sig.start_time + i*sig.time_period] <<"\n";
		out << "[Time:  " << sig.start_time + i * sig.time_period << ", Value: " << sig[sig.start_time + i * sig.time_period] << " ] \n";

	}
	return out;
}


class DS_Filter {
public:
	DS_Filter() {};
	static Signal Low_Pass(Signal x, double dt, double RC) {

		Signal y(x);
		double alpha = dt / (RC + dt);
		y[x.start_time] = alpha * x[x.start_time];

		for (int i = x.start_time + x.time_period; i < x.start_time + x.time_period*x.number_of_samples; i += x.time_period) {
			y[i] = alpha * x[i] + (1 - alpha) * y[i - x.time_period];
		}
		return y;

	}

	static Signal High_Pass(Signal x, double dt, double RC) {

		Signal y(x);
		double alpha = dt / (RC + dt);
		y[x.start_time] = x[x.start_time];



		for (int i = x.start_time+x.time_period; i < x.start_time + x.time_period*x.number_of_samples; i += x.time_period) {
			y[i] = alpha * y[i - x.time_period] + alpha * (x[i] - x[i - x.time_period]);
		}
		return y;

	}
	static Signal Median(Signal x, int window_size) {

		Signal y(x);

		int HWS = std::floor(window_size / 2);



		for (int i = y.start_time; i < y.start_time + y.time_period*HWS; i+=y.time_period) {
			y[i] = x[i];
		}

		for (int i = x.start_time + HWS*x.time_period; i < x.start_time + x.time_period*x.number_of_samples - x.time_period*HWS; i+=x.time_period) {
			std::vector<double> med_sample;
			for (int j = -HWS*x.time_period; j < HWS*x.time_period; j+= x.time_period) {
				med_sample.push_back(x[i + j]);


			}
			double median = Vector_Median(med_sample);
			y[i] = median;


		}
		for (int i = x.start_time + x.number_of_samples*x.time_period - HWS* x.time_period; i < x.start_time + x.time_period*x.number_of_samples; i+= x.time_period) {
			y[i] = x[i];
		}


		return y;

	}
	static Signal Mean(Signal x, int window_size) {

		Signal y(x);

		int HWS = std::floor(window_size / 2);



		for (int i = y.start_time; i < y.start_time + y.time_period*HWS; i += y.time_period) {
			y[i] = x[i];
		}

		for (int i = x.start_time + HWS * x.time_period; i < x.start_time + x.time_period*x.number_of_samples - x.time_period*HWS; i += x.time_period) {
			std::vector<double> med_sample;
			for (int j = -HWS * x.time_period; j < HWS*x.time_period; j += x.time_period) {
				med_sample.push_back(x[i + j]);


			}
			double mean = Vector_Mean(med_sample);
			y[i] = mean;


		}
		for (int i = x.start_time + x.number_of_samples*x.time_period - HWS * x.time_period; i < x.start_time + x.time_period*x.number_of_samples; i += x.time_period) {
			y[i] = x[i];
		}


		return y;

	}

};



class Delta_Signal{
public:
	std::vector<IntDoublePair> Sig;
	int time_period;
	int start_time;
	int loaction_at_time;

	Delta_Signal(int time_period){

		this->time_period = time_period;
		Sig.push_back(IntDoublePair(0, 1));
		loaction_at_time = 0;
		
	};
	Delta_Signal(Delta_Signal const &copy) {
		this->Sig = copy.Sig;
		this->loaction_at_time = copy.loaction_at_time;
		this->start_time = copy.start_time;
		this->time_period = copy.time_period;
	}
	Delta_Signal Flip() {
		Delta_Signal fliped(*this);
		fliped.Sig[0].first *= -1;
		fliped.loaction_at_time = Sig[0].first;
		fliped.start_time = fliped.loaction_at_time;
		return fliped;
	}
	Delta_Signal Shift(int const &units) {
		Delta_Signal Shifted(*this);
		Shifted.Sig[0].first += (-units) * this->time_period;
		Shifted.loaction_at_time = Shifted.Sig[0].first;
		Shifted.start_time = Shifted.Sig[0].first;
		return Shifted;
	}
	double &Value() {
		return this->Sig[0].second;
	}

};
Delta_Signal operator*(double const &value , Delta_Signal const &sig) {
	Delta_Signal result(sig);
	result.Sig[0].second *= value;
	return result;
}
std::ostream &operator<<(std::ostream &out, Delta_Signal const &ds) {
	out << "[Time:  " << ds.Sig[0].first << ", Value: " << ds.Sig[0].second << " ] \n";
	return out;
}

Signal operator+(Delta_Signal &d1, Delta_Signal &d2) {
	assert(d1.time_period == d2.time_period);

	Signal sum(0,std::min(d1.start_time,d2.start_time),d1.time_period);
	if (d1.loaction_at_time < d2.loaction_at_time) {
		sum(d1.Value(), d1.Sig[0].first);
		sum(d2.Value(), d2.Sig[0].first);

		for (int i = d1.loaction_at_time+d1.time_period; i < d2.loaction_at_time; i += d1.time_period) {
			sum(0, i);
		}
	}
	else if (d1.loaction_at_time == d2.loaction_at_time) {
		sum(d1.Value() + d2.Value(), d1.Sig[0].first);
	}
	else {
		sum(d2.Value(), d2.Sig[0].first);
		sum(d1.Value(), d1.Sig[0].first);
		for (int i = d2.loaction_at_time+d1.time_period; i < d1.loaction_at_time; i += d1.time_period) {
			sum(0, i);
		}
	}

	

	
	return sum;
}
Signal operator+(Signal &s, Delta_Signal &ds) {
	assert(s.time_period == ds.time_period);
	Signal res(s);
	try {
		res[0];
		res[0] += ds.Sig[0].second;
	}
	catch (...) {
		res(ds.Sig[0].second, ds.Sig[0].first);
		return res;
	}
	return res;


}


class System {
private:
	std::function<double (double)> func;
public:
	System(std::function<double (double)> expression) {
		this->func = expression;
	}
	double operator()(double const &value) {
		return func(value);
	}
};

Delta_Signal operator*(Delta_Signal &ds, System &sys) {
	Delta_Signal res(ds);
	res.Value() = sys(res.Value());
	return res;
}

Signal operator*(Signal &s, System &sys) {
	Signal res(s);

	for (int i = res.start_time; i < res.start_time + res.number_of_samples*res.time_period; i += res.time_period) {
		res[i] = sys(res[i]);
	}
	return res;
}


std::vector<Complex> Discrete_Fourier_Transform(Signal &dat) {
	std::vector<Complex> Res;
	Res.reserve(dat.number_of_samples);
	int d_size = dat.number_of_samples;


	for (int k = dat.start_time; k < dat.start_time + dat.time_period*dat.number_of_samples; k+=dat.time_period) {
		Complex sum(0, 0);
		for (int n = dat.start_time; n < dat.start_time + dat.time_period*dat.number_of_samples; n += dat.time_period) {
			Complex exp_power(0, ((-2 * std::_Pi) / d_size) *(k*n));
			exp_power = exp_power.Exp();
			exp_power.Real *= dat[n];
			sum = sum + (exp_power);
		}
		Res.emplace_back(sum);
	}
	return Res;
}
Signal Reverse_Discrete_Fourier_Transform(std::vector<Complex> &dat,int const &time_period,int const &start_time) {


	Signal Res(dat.size(), start_time, time_period,DSP_EMPTY);

	int d_size = (int)dat.size();
	for (int k = start_time; k < start_time + time_period* dat.size(); k+=time_period) {
		Complex sum(0, 0);
		for (int n = start_time;n < start_time + time_period * dat.size(); n += time_period) {
			Complex exp_power(0, ((2 * std::_Pi)*(k*n)) / d_size);
			exp_power = exp_power.Exp();
			sum = sum + (dat[n] * exp_power);
		}
		sum = sum * (Complex(1, 0) / Complex(d_size, 0));
		Res[k] = sum.Real;
	}
	return Res;
}






struct Wav_Header {
	char                RIFF[4];        // RIFF Header 
	unsigned long       ChunkSize;      // RIFF Chunk Size  
	char                WAVE[4];        // WAVE Header      
	char                fmt[4];         // FMT header       
	unsigned long       Subchunk1Size;  // Size of the fmt chunk                                
	unsigned short      AudioFormat;    // Audio format 1=PCM,6=mulaw,7=alaw, 257=IBM Mu-Law, 258=IBM A-Law, 259=ADPCM 
	unsigned short      NumOfChan;      // Number of channels 1=Mono 2=Sterio                   
	unsigned long       SamplesPerSec;  // Sampling Frequency in Hz                             
	unsigned long       bytesPerSec;    // bytes per second 
	unsigned short      blockAlign;     // 2=16-bit mono, 4=16-bit stereo 
	unsigned short      bitsPerSample;  // Number of bits per sample      
	char                Subchunk2ID[4]; // "data"  string   
	unsigned long       Subchunk2Size;  // Sampled data length    
};

std::ostream &operator<<(std::ostream &out, Wav_Header const &ds) {
	out << "\n=======================\n";
	out << "RIFF - ";
	for (int i = 0; i < 4; i++)
	{
		out << ds.RIFF[i]<<" ";
	}
	out << "\n";
	out << "RIFF Chunk Size : [ " << ds.ChunkSize << " ]\n";
	out << "WAVE Header - ";
	for (int i = 0; i < 4; i++)
	{
		out << (int)ds.WAVE[i] << " ";
	}
	out << "\nFMT header - ";
	for (int i = 0; i < 4; i++)
	{
		out << (int)ds.fmt[i] << " ";
	}
	out << "\n";
	out << "Size of the fmt chunk : [ " << ds.Subchunk1Size << " ]\n";
	out << "Audio format : [ " << ds.AudioFormat << " ]\n";
	out << "Number of channels : [ " << ds.NumOfChan << " ]\n";
	out << "Sampling Frequency in Hz : [ " << ds.SamplesPerSec << " ]\n";
	out << "bytes per second  : [ " << ds.bytesPerSec << " ]\n";
	out << "2=16-bit mono, 4=16-bit stereo   : [ " << ds.blockAlign << " ]\n";
	out << "Number of bits per sample   : [ " << ds.bitsPerSample << " ]\n";
	out << "data string - ";
	for (int i = 0; i < 4; i++)
	{
		out<< (int)ds.Subchunk2ID[i] <<" ";
	}
	out << "\n";
	out << "Sampled data length  : [ " << ds.Subchunk2Size << " ]\n";


	return out;
}

template <class Word>
std::ostream& write_word(std::ostream& outs, Word value, unsigned size = sizeof(Word))
{
	for (; size; --size, value >>= 8)
		outs.put(static_cast <char> (value & 0xFF));
	return outs;
}


class Wav_File {
public:
	//members
	Wav_Header w_Header;
	unsigned long long f_Size;
	char *Audio_Data;

	//methods
	Wav_File(std::string const &f_path,int mode);
	~Wav_File();
	void Write (std::string const &f_path);
	double get_Duration();
	void add_Sample(std::initializer_list<int> lst);


private:
	std::fstream wavFile;
	int file_mode;
	size_t data_chunk_pos;

};
Wav_File::Wav_File(std::string const &f_path,int mode =0) {
	switch (mode)
	{
	case(0):
		wavFile.open(f_path);
		assert(wavFile.is_open() == true);
		file_mode = mode;
		wavFile.read((char*)&w_Header, sizeof(w_Header));
		unsigned long long begin, end;
		begin = wavFile.tellg();
		wavFile.seekg(0, std::ios::end);
		end = wavFile.tellg();
		f_Size = end - begin;

		wavFile.seekg(44);
		this->Audio_Data = new char[w_Header.Subchunk2Size];
		wavFile.read(Audio_Data, w_Header.Subchunk2Size);
		wavFile.close();
		break;
	case(1):
		wavFile.open(f_path, std::fstream::out |std::ios::binary);
		file_mode = mode;
		wavFile << "RIFF----WAVEfmt ";     // (chunk size to be filled in later)
		write_word(wavFile, 16, 4);  // no extension data
		write_word(wavFile, 1, 2);  // PCM - integer samples
		write_word(wavFile, 1, 2);  // two channels (stereo file)
		write_word(wavFile, 44100, 4);  // samples per second (Hz)
		write_word(wavFile, 176400/2, 4);  // (Sample Rate * BitsPerSample * Channels) / 8
		write_word(wavFile, 4, 2);  // data block size (size of two integer samples, one for each channel, in bytes)
		write_word(wavFile, 16, 2);  // number of bits per sample (use a multiple of 8)
		// Write the data chunk header
		data_chunk_pos = wavFile.tellp();
		wavFile << "data----";  // (chunk size to be filled in later)


		break;
	default:
		break;
	}



}
Wav_File::~Wav_File() {
	if (sizeof(this->Audio_Data) > 0) {
		delete[] this->Audio_Data;
	}
	if (wavFile.is_open() == true) {
		wavFile.close();
	}
}

void Wav_File::Write(std::string const &f_path)
{
	size_t file_length = wavFile.tellp();
	switch (file_mode)
	{
	case(0):
		this->wavFile.open(f_path, std::fstream::out | std::ios::binary);

		wavFile.write((char*)&this->w_Header, sizeof(w_Header));

		wavFile.write(Audio_Data, w_Header.Subchunk2Size);

		this->wavFile.close();
		break;

	case(1):
		// (We'll need the final file size to fix the chunk sizes above)

		// Fix the data chunk header to contain the data size
		wavFile.seekp(data_chunk_pos + 4);
		write_word(wavFile, file_length - data_chunk_pos + 8);

		// Fix the file header to contain the proper RIFF chunk size, which is (file size - 8) bytes
		wavFile.seekp(0 + 4);
		write_word(wavFile, file_length - 8, 4);

		break;
	default:
		break;
	}

}

inline double Wav_File::get_Duration()
{
	return this->w_Header.Subchunk2Size / w_Header.bytesPerSec;
}
void Wav_File::add_Sample(std::initializer_list<int> lst) {
	assert(file_mode == 1);
	for (auto& i : lst)
	{
		write_word(this->wavFile, i, 1);

	}

}
