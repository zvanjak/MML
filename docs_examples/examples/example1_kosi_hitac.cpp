#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Derivation.h"
#endif


using namespace MML;

// TODO 0.9
/*
0. odrediti tocne koordinate!
1. lokalna ravnina za prva tri slucaja, i onda vec za treci Coriolis
2. 3 slucaja gdje je Coriolis presudan
3. 2 slucaja di smo u orbiti
4. 3 slucaja di je specijalna relativnost presudna
*/

void Example1_kosi_hitac()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                     EXAMPLE 1 - kosi hitac                    ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    // lokalni koord sustav i mapiranje na zemljinu kuglu - trazi se lat, long, height
    // 10, 100, 500, 1000, 5000, 10.000, 

    // postaviti lokalni (sferni) koord sustav, fiksiran prema Suncu, u kojem se zemlja vrti, 
    // ali se ovisno o t, tocno zna transf. tocke iz Zemlj. koord u taj sferni sustav

    // imamo pocetnu tocku, za dani T tocno znamu njene sferne koordinate
    
    // Solarni sustav
    // 20.000, 1e5, 1e6, 1e7, 1e8, 2.5e8 m/s, where it will be in 1 hour
}