#ifndef ANALYSIS_PARAMETERS_H
#define ANALYSIS_PARAMETERS_H

#include "Helpers.cc"

// struct Analyzer::Parameters {
// 	OptionalRange<double> pT_small;
// 	OptionalRange<double> pT_large;

// 	OptionalRange<double> y_small;
// 	OptionalRange<double> y_large;

// 	/// Filename of the exported file. 
// 	/// If set to `nullopt`, no data is exported.
// 	std::optional<std::string> filename;

// 	double sigma_eff;
// 	double m;

// 	bool parameter_validation() {
// 		return OptionalRange<double>::disjoint(pT_small, pT_large) || OptionalRange<double>::disjoint(y_small, y_large);
// 	}

// 	Analyzer::Parameters() {
// 		pT_small = OptionalRange<double>();
// 		pT_large = OptionalRange<double>();
// 		y_small = OptionalRange<double>();
// 		y_large = OptionalRange<double>();
// 		filename = std::nullopt;
// 		sigma_eff = Defaults::sigma_eff;
// 		m = Defaults::m;
// 	}

// 	Analyzer::Parameters(
// 		OptionalRange<double> pT_1, 
// 		OptionalRange<double> pT_2, 
// 		OptionalRange<double> y_1, 
// 		OptionalRange<double> y_2, 
// 		std::optional<std::string> fn = std::nullopt,
// 		double _m = Defaults::m,
// 		double seff = Defaults::sigma_eff) {
// 		pT_small = pT_1;
// 		pT_large = pT_2;
// 		y_small = y_1;
// 		y_large = y_2;
// 		filename = fn;
// 		m = _m;
// 		sigma_eff = seff;
// 		assert(parameter_validation());
// 	}

// 	Analyzer::Parameters(
// 		std::optional<double> pT_1_l,
// 		std::optional<double> pT_1_u,
// 		std::optional<double> pT_2_l,
// 		std::optional<double> pT_2_u,
// 		std::optional<double> y_1_l,
// 		std::optional<double> y_1_u,
// 		std::optional<double> y_2_l,
// 		std::optional<double> y_2_u,
// 		std::optional<string> fn = std::nullopt,
// 		double _m = Defaults::m,
// 		double seff = Defaults::sigma_eff) {
// 		pT_small = OptionalRange<double>(pT_1_l, pT_1_u);
// 		pT_large = OptionalRange<double>(pT_2_l, pT_2_u);
// 		y_small = OptionalRange<double>(y_1_l, y_1_u);
// 		y_large = OptionalRange<double>(y_2_l, y_2_u);
// 		filename = fn;
// 		m = _m;
// 		sigma_eff = seff;
// 		assert(parameter_validation());
// 	}
// };

#endif // ANALYSIS_PARAMTERS_H