#include <iostream>
#include <fstream>
#include <vector>

using v_t = std::vector <double>;

namespace data1
{
	double ro_p, T_b, L, m_v, T_inf, P_a, m_a, c_p, X_inf, D_p;
	double lambda_g, Y_inf;
}

namespace data2
{
	double D, P, Q, S, Y_oa, c_p, lambda_g, T_a, T_b, m_b, L;
}

namespace results
{
	double T_w, X_w, Y_w, ro_w, v_w, m_point_slash, m_point, t_0;
	double R_star;
}

auto read_data()
{
	std::ifstream fin("data.txt"); bool flag = false; std::string data;
	std::vector <double> v; v.reserve(21);

	while (fin >> data)
	{
		if (flag) { v.push_back(std::atof(data.c_str())); flag = false; }
		if (data == "=") flag = true;
	};

	data1::ro_p = v[0]; data1::T_b = v[1]; data1::L = v[2]; data1::m_v = v[3];
	data1::T_inf = v[4]; data1::P_a = v[5]; data1::m_a = v[6]; data1::c_p = v[7];
	data1::X_inf = v[8]; data1::D_p = v[9];

	data2::D = v[10] * 1e-3; data2::P = v[11]; data2::Q = v[12] * 1e3;
	data2::S = v[13]; data2::Y_oa = v[14]; data2::c_p = v[15];
	data2::lambda_g = v[16]; data2::T_a = v[17]; data2::T_b = v[18];
	data2::m_b = v[19] * 1e-3; data2::L = v[20];
}

auto calculate_results_1()
{
	const double R = 8.314462;
	const double kb = 1.380649 * 1e-23;
	const double Na = 6.022 * 1e23;
	const double pi = 3.1416;

	data1::Y_inf = data1::X_inf * data1::m_v / (data1::m_v * data1::X_inf + (1 - data1::X_inf) * data1::m_a);
	{
		double eps = 1; double t = data1::T_b;
		while (true)
		{
			results::X_w = std::exp(data1::L * data1::m_v * 1e-3 / R * (1 / data1::T_b - 1 / t));
			results::Y_w = results::X_w * data1::m_v / (data1::m_v * results::X_w + (1 - results::X_w) * data1::m_a);
			results::T_w = data1::T_inf - data1::L * (results::Y_w - data1::Y_inf) / (data1::c_p * (1 - results::Y_w));

			if (std::abs(results::T_w - t) < eps) { results::T_w = t; break; }
			else t -= 0.001;
		}
	}

	double mu = (results::X_w * data1::m_v + (1. - results::X_w) * data1::m_a) * 1e-3;
	results::ro_w = data1::P_a * mu / (results::T_w * R);
	double v = std::sqrt(3. * R * results::T_w / mu);
	double ro = data1::P_a * mu / (results::T_w * R);
	double d = 2*std::cbrt(3. * mu / (Na * 4. * pi * ro));
	double l = kb * results::T_w / (std::sqrt(2.) * pi * data1::P_a * d * d);
	data1::lambda_g = 0.04635537370193092;

	results::m_point_slash = data1::lambda_g / (data1::c_p * data1::D_p * 1e-3 / 2) * std::log(1. + data1::c_p * (data1::T_inf - results::T_w) / data1::L) * 0.1;
	results::v_w = results::m_point_slash * 1e4 / (results::ro_w);
	results::m_point = results::m_point_slash * pi * data1::D_p * data1::D_p * 1e-2;
	results::t_0 = data1::ro_p * data1::c_p * (data1::D_p * 1e-3 / 2) * (data1::D_p * 1e-3 / 2) / (2. * data1::lambda_g * std::log(1. + data1::c_p * (data1::T_inf - results::T_w) / data1::L));
}

auto write_results()
{
	std::remove("result.txt"); std::ofstream fout("result.txt");
	fout << "ЗАДАНИЕ 1" << std::endl;
	fout << "Темпертура поверхности T_w: " << results::T_w << " К" << std::endl;
	fout << "Мольная доля пара на поверхности X_w: " << results::X_w << std::endl;
	fout << "Массовая доля пара на поверхности Y_w: " << results::Y_w << std::endl;
	fout << "Плотность газа на поверхности ro_w: " << results::ro_w << " кг/м^3" << std::endl;
	fout << "Скорость пара на поверхности v_w: " << results::v_w << " мм/c" << std::endl;
	fout << "Удельный массовый поток с единицы поверхности m`\": " << results::m_point_slash << " г/(см^2*с)" << std::endl;
	fout << "Массовая скорость испарения m` в начальный момент времени: " << results::m_point << " г/с" << std::endl;
	fout << "Время полного испарения капли t_o: " << results::t_0 << " с" << std::endl;
	fout << "Tеплопроводность газа lambda_g: " << data1::lambda_g << " Вт/(м*K)" << std::endl << std::endl;
	fout << "ЗАДАНИЕ 2" << std::endl;
	fout << "Радиус фронта пламени R_star: " << results::R_star << " мм";
}

auto write_radial_distributions_1(v_t& Yv_r, v_t& T_r, v_t& vg_r, v_t& r)
{
	std::remove("radial_distributions_1.txt");
	std::ofstream fout("radial_distributions_1.txt");
	fout << "r, мм	-	Yv		-	Tg, К	-	vg, мм/с" << std::endl;
	for (auto i = 0U; i < std::size(r); i++)
		fout << r[i] << "		" << Yv_r[i] << "   		" << T_r[i] << "		" << vg_r[i] << std::endl;
}

auto radial_distributions_1()
{
	const double R = 8.314462;
	std::size_t n = 1e4;
	double dr = 19 * data1::D_p / (2. * n);
	v_t Yv_r; v_t T_r; v_t vg_r; v_t r;
	Yv_r.reserve(n); T_r.reserve(n); vg_r.reserve(n); r.reserve(n);

	double b = data1::c_p * results::ro_w * results::v_w * data1::D_p * data1::D_p * 1e-6 / (4. * data1::lambda_g);
	double point = data1::D_p / 2.;
	double c1 = (data1::T_inf - results::T_w) / (1. - std::exp(-b / point));
	double c2 = data1::T_inf - c1;

	double a = results::ro_w * results::v_w * data1::D_p * data1::D_p / 4.;

	for (auto i = 0U; i < n; i++)
	{
		r.push_back(point);
		T_r.push_back(c1 * std::exp(-b / point) + c2);

		double t = data1::c_p * (data1::T_inf - T_r.back()) / data1::L;
		Yv_r.push_back((t + data1::Y_inf) / (1. + t));

		double X = Yv_r.back() / data1::m_v / (Yv_r.back() / data1::m_v + (1 - Yv_r.back()) / data1::m_a);

		double mu = (X * data1::m_v + (1. - X) * data1::m_a) * 1e-3;
		double ro = data1::P_a * mu / (T_r.back() * R);
		vg_r.push_back(a / (ro * point * point));

		point += dr;
	}

	write_radial_distributions_1(Yv_r, T_r, vg_r, r);
}

auto write_radial_distributions_2(v_t& Yf_r, v_t& Yo_r, v_t& T_r, v_t& r)
{
	std::remove("radial_distributions_2.txt");
	std::ofstream fout("radial_distributions_2.txt");
	fout << "r, мм	-	Yf_r		-	Yo_r	-	Tr, К" << std::endl;
	for (auto i = 0U; i < std::size(r); i++)
		fout << r[i] << "		" << Yf_r[i] << "   		" << Yo_r[i] << "		" << T_r[i] << std::endl;
}

auto calculate_results_2()
{
	const double P_a = 101325;
	const double R = 8.314462;
	const double pi = 3.1416;
	double T_w = 1 / (1 / data2::T_b - std::log(data2::P / P_a) * data2::m_b * R / data2::L);
	double t = (data2::c_p * (data2::T_a - T_w) + data2::Q * data2::Y_oa / data2::S) / data2::L;
	double m_point = 2 * pi * data2::D * data2::lambda_g * std::log(1 + t) / data2::c_p;
	double F_w = 1 - (1 + data2::Y_oa / data2::lambda_g) / (1 + t);
	double b = m_point * data2::c_p * data2::D * data2::D / (4 * data2::lambda_g);
	double c1 = -1 / (1 - std::exp(-b * 2 / data2::D));
	double c2 = -c1;
	b *= 1e3;

	results::R_star = -b / std::log((data2::Y_oa / (data2::S * F_w + data2::Y_oa) - c2) / c1);

	std::size_t n = 1e4;
	double point = data2::D / 2 * 1e3;
	double max_point = point + (results::R_star - point) * 3;
	double dr = (max_point - point) / n;
	v_t Yf_r; v_t Yo_r; v_t T_r; v_t r;
	Yf_r.reserve(n); Yo_r.reserve(n); T_r.reserve(n); r.reserve(n);

	double beta_w = F_w;
	double beta_a = -data2::Y_oa / data2::S;
	double beta_tw = data2::c_p * T_w + data2::Q * F_w;
	double beta_ta = data2::c_p * data2::T_a;

	while (point <= max_point)
	{
		r.push_back(point);
		double beta_tild = c1 * std::exp(-b / r.back()) + c2;
		double beta = beta_tild * (beta_w - beta_a) + beta_a;
		double beta_t = beta_tild * (beta_tw - beta_ta) + beta_ta;

		if (point <= results::R_star)
		{
			Yo_r.push_back(0);
			Yf_r.push_back(beta);
			T_r.push_back((beta_t - data2::Q * Yf_r.back()) / data2::c_p);
		}
		else
		{
			Yf_r.push_back(0);
			Yo_r.push_back(-data2::S * beta);
			T_r.push_back(beta_t / data2::c_p);
		}
		point += dr;
	}

	write_radial_distributions_2(Yf_r, Yo_r, T_r, r);
}

int main()
{
	read_data();
	calculate_results_1();
	radial_distributions_1();
	calculate_results_2();
	write_results();
}