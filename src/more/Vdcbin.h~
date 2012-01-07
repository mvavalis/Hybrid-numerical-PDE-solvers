#ifndef _VDCBIN_H_
#define _VDCBIN_H_

/* set MSB to 1, all other bits to 0 */
#define ULLONG_MSB_1 (~((((unsigned long long)-1)<<1)>>1))
#ifndef ULLONG_MAX
	#define ULLONG_MAX (~((unsigned long long)0))
#endif

/**
 * Gray code van der Corput sequence:
 * How? (when the word "bit" is used it may mean the bit's position)
 * v: van der Corput sequence
 * d: driver
 * constant MSB_1: MSB is set to 1, all other bits are set to 0
 *
 * STEP 0:
 *   [v=0; d=0;]
 * STEP 1: increment d and find the bit that was altered from 0 to 1 (that's only one bit)
 *   [tmp = ~d & ++d;]
 * STEP 2: mirror this bit (e.g.: for int length equal to 7: 0000010 -> 0100000)
 *   [tmp = MSB_1/tmp;]
 * STEP 3: alter the mirrored bit in v
 *   [v ^= tmp;]
 */
class Vdcbin
{
private:
	unsigned long long v, d; //van der Corput sequence, driver
	
	inline void next() {v ^= ULLONG_MSB_1/(~d & ++d);}
	inline void next_cz() { unsigned long long tmp; //check div by 0
		v ^= ULLONG_MSB_1/(((tmp = ~d & ++d) == 0) ? d=1:tmp);}
public:
	Vdcbin(unsigned long long v=0) : d(0) {this->v=v;}
	inline unsigned long long ullint()    {next();    return v;}
	inline unsigned long long ullint_cz() {next_cz(); return v;}
	inline double doub()    {next();    return v/(ULLONG_MAX+1.);}
	inline double doub_cz() {next_cz(); return v/(ULLONG_MAX+1.);}
};

#endif
