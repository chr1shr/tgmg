#ifndef TGMGPP_PREDICT_HH
#define TGMGPP_PREDICT_HH

#include "tgmg_config.hh"

struct tgmg_predict {
    /** A multiplier to apply to the internal counter for tracking the number
     * of V-cycles required. */
	const int mult;
    /** A decay factor to try and lower the number of multigrid cycles
     * performed at each time. */
	const int decay;
    /** The threshold on the value of lim, which sets a limit on the number
     * of V-cycles performed before bailing out with an error. */
	const int max_thresh;
    /** The number of extra iterations to perform after reaching the required
     * tolerance. */
	const int extra_iters;
	/* A current estimate of the number of V-cycles required for
	 * convergence, multiplied by the mgs_mult constant. */
	int lim;
    /** A running counter of the number of multigrid solves that have been
     * performed, used for printing diagnostic information about the average
     * number of V-cycles. */
    int solves;
	/** A running counter of the number of multigrid V-cycles that have
	 * been performed. This does not count the number of extra iterations. */
	int vcount;
    /** Initializes the tgmg_predict class that is used to predict the number
     * of V-cycles required to reach the tolerance level.
     * \param[in] extra_iters the number of extra iterations to apply after
     *                        reaching the tolerance level. */
	tgmg_predict(int extra_iters_=tgmg_predict_extra_iters)
        : mult(tgmg_predict_mult), decay(tgmg_predict_decay),
		max_thresh(mult*tgmg_predict_max_iters), extra_iters(extra_iters_),
		lim(mult*tgmg_predict_init), solves(0), vcount(0) {}
    /** Initializes the tgmg_predict class that is used to predict the number
     * of V-cycles required to reach the tolerance level.
     * \param[in] (mult_,decay_,max_iters,init_iters) the prediction parameters.
     * \param[in] extra_iters the number of extra iterations to apply after
     *                        reaching the tolerance level. */
	tgmg_predict(int mult_,int decay_,int max_iters,int init_iters,int extra_iters_=tgmg_predict_extra_iters)
        : mult(mult_), decay(decay_), max_thresh(mult*max_iters),
        extra_iters(extra_iters_), lim(mult*init_iters), solves(0),
        vcount(0) {}
    /** Resets the predicted number of V-cycles. */
    inline void reset() {
        lim=mult*tgmg_predict_init;
    }
    /** Updates the counters for the number of solves, and the number V-cycles
     * that have been performed.
     * \param[in] n the number of V-cycles (not counting extra iterations) that
     *              were performed on the current solve. */
    inline void add_iters(int n) {
		solves++;
		vcount+=n;
	}
    /** Returns the total number of V-cycles performed, including extra
     * iterations.
     * \return The total number. */
    inline int vcount_extra() {
        return vcount+extra_iters*solves;
    }
    /** Returns the average number of V-cycles based on internal counters in
     * the tgmg_predict class. Those counters are reset in this call.
     * \param[in] count_extra whether to count the extra iterations as well.
     * \return The average number of V-cycles. */
	inline float avg_iters(bool count_extra=false) {
		float ans=static_cast<float>(vcount)/static_cast<float>(solves);
		solves=vcount=0;
		return count_extra?ans+extra_iters:ans;
	}
};

#endif
