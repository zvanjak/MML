//  __ __ _____ _
// |  |  |  __ | |    Minimal Physics Library for Modern C++
// | | | |  ___| |__  version 0.6
// |_| |_|_|   |____| https://github.com/zvanjak/mpl
//
// Copyright: 2025 - , Zvonimir Vanjak 
//
// LICENSING:
// - GPL 3.0 license for personal, educational and research use.
// - Commercial license available upon request

#if !defined MPL_BASE_H
#define MPL_BASE_H


namespace MPL
{
	enum CollisionSimulatorRunType
	{
		RunTypeExact = 0,					// exact simulation, 
		RunTypeFast,							// fast simulation, with subdivision of the box into subcontainers,
		RunTypeFastMultithread,		// fast simulation, with subdivision and using multithreading
		RunTypeFastMultithreadTP	// fast simulation, with subdivision and using multithreading with thread pool
	};
}

#endif // MPLBase.h