#pragma once

class GaussianDistribution {
public:
	GaussianDistribution() {}

	double getMu() const {
		return mu;
	}

	double getSigma() {
		if(needUpdate) {
			sigma = sqrt(sn / datapoints);
			needUpdate = false;
		}
		return sigma;
	}

	void addDataPoint(double value) {
		datapoints++;
		double oldMu = mu;
		mu += (value - mu) / datapoints;
		sn += (value - oldMu) * (value - mu);
		needUpdate = true;
		assert(!std::isnan(mu) && !std::isinf(mu) && !std::isnan(sn) && !std::isinf(sn));
	}

	void removeDataPoint(double value) {
		datapoints--;
		if(datapoints == 0) {
			mu = 0;
			sn = 0;
			sigma = 0;
			return;
		}
		double oldMu = mu;
		mu -= (value - mu) / datapoints;
		sn -= (value - oldMu) * (value - mu);
		needUpdate = true;
		assert(!std::isnan(mu) && !std::isinf(mu) && !std::isnan(sn) && !std::isinf(sn));
	}

	double getCDF(double value) {
		return 0.5 * (1 + std::erf((value - mu) / (sqrt(2 * getSigma()))));
	}

	GaussianDistribution operator+(GaussianDistribution a) {
		GaussianDistribution gd;
		gd.mu = mu + a.mu;
		gd.sigma = getSigma() + getSigma();
		return gd;
	}

	GaussianDistribution operator*(double v) {
		GaussianDistribution gd;
		gd.mu = v * mu;
		gd.sigma = getSigma();
		return gd;
	}

protected:
	double mu = 0, sn = 0, sigma = 0;
	double datapoints = 0;
	bool needUpdate = false;
};