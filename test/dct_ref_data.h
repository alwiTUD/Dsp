#pragma once
#include <vector>
// Define dct test input - output pairs for all allowed data types.
	// signal length = n
std::vector<long double> xLongDouble{ 4.0, 2.0, 8.0, 9.0, 4.0, 7.0, 3.0, 6.0, 9.0, 6.0 };

std::vector<double> xDouble{ 4.0, 2.0, 8.0, 9.0, 4.0, 7.0, 3.0, 6.0, 9.0, 6.0 };

std::vector<float> xFloat{ 4.0f, 2.0f, 8.0f, 9.0f, 4.0f, 7.0f, 3.0f, 6.0f, 9.0f, 6.0f };

// signal length = n
const std::vector<long double> yLongDouble{ 18.341210428976602,
										  -2.031945911545793,
										  -0.688190960235587,
										  -3.960267084672950,
										  -1.841640786499874,
										   1.897366596101028,
										   0.162459848116453,
										   4.933516074140680,
										   0.841640786499874,
										  -1.802192990998431 };


const std::vector<double> yDouble{ 18.341210428976602,
								  -2.031945911545793,
								  -0.688190960235587,
								  -3.960267084672950,
								  -1.841640786499874,
								   1.897366596101028,
								   0.162459848116453,
								   4.933516074140680,
								   0.841640786499874,
								  -1.802192990998431 };

const std::vector<float> yFloat{ 18.341210428976602f,
							  -2.031945911545793f,
							  -0.688190960235587f,
							  -3.960267084672950f,
							  -1.841640786499874f,
							   1.897366596101028f,
							   0.162459848116453f,
							   4.933516074140680f,
							   0.841640786499874f,
							  -1.802192990998431f };

// n = 5 (less than signal length, dim. reduction)
std::vector<long double> yLongDouble5 = { 12.074767078498866,
										  -2.602236241221291,
										  -3.116140650889268,
										   4.210506685052820,
										   0.994820307329625 };

std::vector<long double> yDouble5 = { 12.074767078498866,
									  -2.602236241221291,
									  -3.116140650889268,
									   4.210506685052820,
									   0.994820307329625 };

std::vector<long double> yFloat5 = { 12.074767078498866f,
									  -2.602236241221291f,
									  -3.116140650889268f,
									   4.210506685052820f,
									   0.994820307329625f };

// n = 15 (more than signal length, added dimensions).
std::vector<long double> yLongDouble15 = { 14.975535605335345,
										   7.670713755759988,
										  -5.997901881420042,
										  -0.561905566057725,
										   0.329917647780843,
										  -5.692099788303081,
										  -1.503693405474195,
										   2.955102356260902,
										  -0.530850817028425,
										   0.132647910525122,
										   4.016632088371217,
										   2.878859482290704,
										   0.687196824546469,
										  -0.637460029419356,
										  -1.871147092899867 };

std::vector<long double> yDouble15 = { 14.975535605335345,
									   7.670713755759988,
									  -5.997901881420042,
									  -0.561905566057725,
									   0.329917647780843,
									  -5.692099788303081,
									  -1.503693405474195,
									   2.955102356260902,
									  -0.530850817028425,
									   0.132647910525122,
									   4.016632088371217,
									   2.878859482290704,
									   0.687196824546469,
									  -0.637460029419356,
									  -1.871147092899867 };

std::vector<long double> yFloat15 = { 14.975535605335345f,
								   7.670713755759988f,
								  -5.997901881420042f,
								  -0.561905566057725f,
								   0.329917647780843f,
								  -5.692099788303081f,
								  -1.503693405474195f,
								   2.955102356260902f,
								  -0.530850817028425f,
								   0.132647910525122f,
								   4.016632088371217f,
								   2.878859482290704f,
								   0.687196824546469f,
								  -0.637460029419356f,
								  -1.871147092899867f };

// Define test target(s) for the cosineBasisVector calculation.
const std::vector<std::vector<long double>> basisVectorsLongDouble =
{
	{0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979},
	{0.315252941349890, 0.307490367669328, 0.292156360634725, 0.269628494389924, 0.240461479767574, 0.205373505467449, 0.165228553867129, 0.121015126908468, 0.0738219058991409, 0.0248109445657177, -0.0248109445657176, -0.0738219058991408, -0.121015126908468, -0.165228553867129, -0.205373505467449, -0.240461479767574, -0.269628494389924, -0.292156360634725, -0.307490367669328, -0.315252941349890},
	{0.312334477467278, 0.281761002650515, 0.223606797749979, 0.143564401525505, 0.0494689214077114, -0.0494689214077113, -0.143564401525505, -0.223606797749979, -0.281761002650515, -0.312334477467278, -0.312334477467278, -0.281761002650515, -0.223606797749979, -0.143564401525505, -0.0494689214077114, 0.0494689214077113, 0.143564401525505, 0.223606797749979, 0.281761002650515, 0.312334477467278},
	{0.307490367669328, 0.240461479767574, 0.121015126908468, -0.0248109445657176, -0.165228553867129, -0.269628494389924, -0.315252941349890, -0.292156360634725, -0.205373505467449, -0.0738219058991411, 0.0738219058991409, 0.205373505467449, 0.292156360634725, 0.315252941349890, 0.269628494389924, 0.165228553867129, 0.0248109445657176, -0.121015126908468, -0.240461479767574, -0.307490367669328},
	{0.300750477503773, 0.185874017230092, 1.93633660727019e-17, -0.185874017230092, -0.300750477503773, -0.300750477503773, -0.185874017230092, -5.80900982181058e-17, 0.185874017230092, 0.300750477503773, 0.300750477503773, 0.185874017230092, 9.68168303635097e-17, -0.185874017230092, -0.300750477503773, -0.300750477503773, -0.185874017230092, -1.35543562508914e-16, 0.185874017230092, 0.300750477503773},
	{0.292156360634725, 0.121015126908468, -0.121015126908468, -0.292156360634725, -0.292156360634725, -0.121015126908468, 0.121015126908468, 0.292156360634725, 0.292156360634725, 0.121015126908468, -0.121015126908468, -0.292156360634725, -0.292156360634725, -0.121015126908468, 0.121015126908468, 0.292156360634725, 0.292156360634725, 0.121015126908468, -0.121015126908467, -0.292156360634725},
	{0.281761002650515, 0.0494689214077114, -0.223606797749979, -0.312334477467278, -0.143564401525505, 0.143564401525504, 0.312334477467278, 0.223606797749979, -0.0494689214077115, -0.281761002650515, -0.281761002650515, -0.0494689214077112, 0.223606797749979, 0.312334477467278, 0.143564401525505, -0.143564401525504, -0.312334477467278, -0.223606797749979, 0.0494689214077109, 0.281761002650515},
	{0.269628494389924, -0.0248109445657176, -0.292156360634725, -0.240461479767574, 0.0738219058991409, 0.307490367669328, 0.205373505467449, -0.121015126908468, -0.315252941349890, -0.165228553867129, 0.165228553867129, 0.315252941349890, 0.121015126908468, -0.205373505467449, -0.307490367669328, -0.0738219058991409, 0.240461479767574, 0.292156360634725, 0.0248109445657175, -0.269628494389924},
	{0.255833636800846, -0.0977197537924274, -0.316227766016838, -0.0977197537924275, 0.255833636800846, 0.255833636800847, -0.0977197537924273, -0.316227766016838, -0.0977197537924275, 0.255833636800846, 0.255833636800847, -0.0977197537924272, -0.316227766016838, -0.0977197537924276, 0.255833636800846, 0.255833636800847, -0.0977197537924272, -0.316227766016838, -0.0977197537924277, 0.255833636800846},
	{0.240461479767574, -0.165228553867129, -0.292156360634725, 0.0738219058991409, 0.315252941349890, 0.0248109445657179, -0.307490367669328, -0.121015126908468, 0.269628494389924, 0.205373505467450, -0.205373505467449, -0.269628494389924, 0.121015126908468, 0.307490367669328, -0.0248109445657181, -0.315252941349890, -0.0738219058991410, 0.292156360634724, 0.165228553867129, -0.240461479767574},
	{0.223606797749979, -0.223606797749979, -0.223606797749979, 0.223606797749979, 0.223606797749979, -0.223606797749979, -0.223606797749979, 0.223606797749979, 0.223606797749979, -0.223606797749979, -0.223606797749979, 0.223606797749979, 0.223606797749979, -0.223606797749979, -0.223606797749979, 0.223606797749979, 0.223606797749978, -0.223606797749979, -0.223606797749980, 0.223606797749979},
	{0.205373505467449, -0.269628494389924, -0.121015126908468, 0.307490367669328, 0.0248109445657176, -0.315252941349890, 0.0738219058991411, 0.292156360634725, -0.165228553867129, -0.240461479767575, 0.240461479767574, 0.165228553867129, -0.292156360634725, -0.0738219058991421, 0.315252941349890, -0.0248109445657169, -0.307490367669328, 0.121015126908467, 0.269628494389924, -0.205373505467448},
	{0.185874017230092, -0.300750477503773, -5.80900982181058e-17, 0.300750477503773, -0.185874017230092, -0.185874017230093, 0.300750477503773, 7.36003649626590e-16, -0.300750477503773, 0.185874017230092, 0.185874017230092, -0.300750477503773, -8.52183846062801e-16, 0.300750477503773, -0.185874017230092, -0.185874017230093, 0.300750477503773, -1.55102667445532e-16, -0.300750477503773, 0.185874017230092},
	{0.165228553867129, -0.315252941349890, 0.121015126908468, 0.205373505467450, -0.307490367669328, 0.0738219058991406, 0.240461479767574, -0.292156360634725, 0.0248109445657181, 0.269628494389924, -0.269628494389924, -0.0248109445657175, 0.292156360634725, -0.240461479767574, -0.0738219058991411, 0.307490367669328, -0.205373505467450, -0.121015126908469, 0.315252941349890, -0.165228553867129},
	{0.143564401525505, -0.312334477467278, 0.223606797749979, 0.0494689214077117, -0.281761002650515, 0.281761002650515, -0.0494689214077120, -0.223606797749980, 0.312334477467278, -0.143564401525504, -0.143564401525504, 0.312334477467278, -0.223606797749979, -0.0494689214077125, 0.281761002650516, -0.281761002650515, 0.0494689214077106, 0.223606797749980, -0.312334477467278, 0.143564401525503},
	{0.121015126908468, -0.292156360634725, 0.292156360634725, -0.121015126908468, -0.121015126908468, 0.292156360634725, -0.292156360634725, 0.121015126908468, 0.121015126908468, -0.292156360634725, 0.292156360634725, -0.121015126908467, -0.121015126908469, 0.292156360634725, -0.292156360634724, 0.121015126908467, 0.121015126908468, -0.292156360634725, 0.292156360634725, -0.121015126908468},
	{0.0977197537924274, -0.255833636800846, 0.316227766016838, -0.255833636800846, 0.0977197537924273, 0.0977197537924281, -0.255833636800847, 0.316227766016838, -0.255833636800846, 0.0977197537924271, 0.0977197537924277, -0.255833636800847, 0.316227766016838, -0.255833636800846, 0.0977197537924270, 0.0977197537924279, -0.255833636800847, 0.316227766016838, -0.255833636800846, 0.0977197537924268},
	{0.0738219058991409, -0.205373505467450, 0.292156360634725, -0.315252941349890, 0.269628494389924, -0.165228553867128, 0.0248109445657181, 0.121015126908469, -0.240461479767574, 0.307490367669328, -0.307490367669328, 0.240461479767574, -0.121015126908468, -0.0248109445657177, 0.165228553867129, -0.269628494389925, 0.315252941349890, -0.292156360634724, 0.205373505467449, -0.0738219058991401},
	{0.0494689214077114, -0.143564401525505, 0.223606797749979, -0.281761002650515, 0.312334477467278, -0.312334477467278, 0.281761002650515, -0.223606797749979, 0.143564401525505, -0.0494689214077107, -0.0494689214077114, 0.143564401525505, -0.223606797749979, 0.281761002650515, -0.312334477467278, 0.312334477467278, -0.281761002650515, 0.223606797749978, -0.143564401525503, 0.0494689214077104},
	{0.0248109445657177, -0.0738219058991411, 0.121015126908468, -0.165228553867129, 0.205373505467450, -0.240461479767575, 0.269628494389924, -0.292156360634725, 0.307490367669328, -0.315252941349890, 0.315252941349890, -0.307490367669328, 0.292156360634725, -0.269628494389923, 0.240461479767575, -0.205373505467449, 0.165228553867128, -0.121015126908468, 0.0738219058991390, -0.0248109445657143}
};

const std::vector<std::vector<double>> basisVectorsDouble =
{
	{0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979},
	{0.315252941349890, 0.307490367669328, 0.292156360634725, 0.269628494389924, 0.240461479767574, 0.205373505467449, 0.165228553867129, 0.121015126908468, 0.0738219058991409, 0.0248109445657177, -0.0248109445657176, -0.0738219058991408, -0.121015126908468, -0.165228553867129, -0.205373505467449, -0.240461479767574, -0.269628494389924, -0.292156360634725, -0.307490367669328, -0.315252941349890},
	{0.312334477467278, 0.281761002650515, 0.223606797749979, 0.143564401525505, 0.0494689214077114, -0.0494689214077113, -0.143564401525505, -0.223606797749979, -0.281761002650515, -0.312334477467278, -0.312334477467278, -0.281761002650515, -0.223606797749979, -0.143564401525505, -0.0494689214077114, 0.0494689214077113, 0.143564401525505, 0.223606797749979, 0.281761002650515, 0.312334477467278},
	{0.307490367669328, 0.240461479767574, 0.121015126908468, -0.0248109445657176, -0.165228553867129, -0.269628494389924, -0.315252941349890, -0.292156360634725, -0.205373505467449, -0.0738219058991411, 0.0738219058991409, 0.205373505467449, 0.292156360634725, 0.315252941349890, 0.269628494389924, 0.165228553867129, 0.0248109445657176, -0.121015126908468, -0.240461479767574, -0.307490367669328},
	{0.300750477503773, 0.185874017230092, 1.93633660727019e-17, -0.185874017230092, -0.300750477503773, -0.300750477503773, -0.185874017230092, -5.80900982181058e-17, 0.185874017230092, 0.300750477503773, 0.300750477503773, 0.185874017230092, 9.68168303635097e-17, -0.185874017230092, -0.300750477503773, -0.300750477503773, -0.185874017230092, -1.35543562508914e-16, 0.185874017230092, 0.300750477503773},
	{0.292156360634725, 0.121015126908468, -0.121015126908468, -0.292156360634725, -0.292156360634725, -0.121015126908468, 0.121015126908468, 0.292156360634725, 0.292156360634725, 0.121015126908468, -0.121015126908468, -0.292156360634725, -0.292156360634725, -0.121015126908468, 0.121015126908468, 0.292156360634725, 0.292156360634725, 0.121015126908468, -0.121015126908467, -0.292156360634725},
	{0.281761002650515, 0.0494689214077114, -0.223606797749979, -0.312334477467278, -0.143564401525505, 0.143564401525504, 0.312334477467278, 0.223606797749979, -0.0494689214077115, -0.281761002650515, -0.281761002650515, -0.0494689214077112, 0.223606797749979, 0.312334477467278, 0.143564401525505, -0.143564401525504, -0.312334477467278, -0.223606797749979, 0.0494689214077109, 0.281761002650515},
	{0.269628494389924, -0.0248109445657176, -0.292156360634725, -0.240461479767574, 0.0738219058991409, 0.307490367669328, 0.205373505467449, -0.121015126908468, -0.315252941349890, -0.165228553867129, 0.165228553867129, 0.315252941349890, 0.121015126908468, -0.205373505467449, -0.307490367669328, -0.0738219058991409, 0.240461479767574, 0.292156360634725, 0.0248109445657175, -0.269628494389924},
	{0.255833636800846, -0.0977197537924274, -0.316227766016838, -0.0977197537924275, 0.255833636800846, 0.255833636800847, -0.0977197537924273, -0.316227766016838, -0.0977197537924275, 0.255833636800846, 0.255833636800847, -0.0977197537924272, -0.316227766016838, -0.0977197537924276, 0.255833636800846, 0.255833636800847, -0.0977197537924272, -0.316227766016838, -0.0977197537924277, 0.255833636800846},
	{0.240461479767574, -0.165228553867129, -0.292156360634725, 0.0738219058991409, 0.315252941349890, 0.0248109445657179, -0.307490367669328, -0.121015126908468, 0.269628494389924, 0.205373505467450, -0.205373505467449, -0.269628494389924, 0.121015126908468, 0.307490367669328, -0.0248109445657181, -0.315252941349890, -0.0738219058991410, 0.292156360634724, 0.165228553867129, -0.240461479767574},
	{0.223606797749979, -0.223606797749979, -0.223606797749979, 0.223606797749979, 0.223606797749979, -0.223606797749979, -0.223606797749979, 0.223606797749979, 0.223606797749979, -0.223606797749979, -0.223606797749979, 0.223606797749979, 0.223606797749979, -0.223606797749979, -0.223606797749979, 0.223606797749979, 0.223606797749978, -0.223606797749979, -0.223606797749980, 0.223606797749979},
	{0.205373505467449, -0.269628494389924, -0.121015126908468, 0.307490367669328, 0.0248109445657176, -0.315252941349890, 0.0738219058991411, 0.292156360634725, -0.165228553867129, -0.240461479767575, 0.240461479767574, 0.165228553867129, -0.292156360634725, -0.0738219058991421, 0.315252941349890, -0.0248109445657169, -0.307490367669328, 0.121015126908467, 0.269628494389924, -0.205373505467448},
	{0.185874017230092, -0.300750477503773, -5.80900982181058e-17, 0.300750477503773, -0.185874017230092, -0.185874017230093, 0.300750477503773, 7.36003649626590e-16, -0.300750477503773, 0.185874017230092, 0.185874017230092, -0.300750477503773, -8.52183846062801e-16, 0.300750477503773, -0.185874017230092, -0.185874017230093, 0.300750477503773, -1.55102667445532e-16, -0.300750477503773, 0.185874017230092},
	{0.165228553867129, -0.315252941349890, 0.121015126908468, 0.205373505467450, -0.307490367669328, 0.0738219058991406, 0.240461479767574, -0.292156360634725, 0.0248109445657181, 0.269628494389924, -0.269628494389924, -0.0248109445657175, 0.292156360634725, -0.240461479767574, -0.0738219058991411, 0.307490367669328, -0.205373505467450, -0.121015126908469, 0.315252941349890, -0.165228553867129},
	{0.143564401525505, -0.312334477467278, 0.223606797749979, 0.0494689214077117, -0.281761002650515, 0.281761002650515, -0.0494689214077120, -0.223606797749980, 0.312334477467278, -0.143564401525504, -0.143564401525504, 0.312334477467278, -0.223606797749979, -0.0494689214077125, 0.281761002650516, -0.281761002650515, 0.0494689214077106, 0.223606797749980, -0.312334477467278, 0.143564401525503},
	{0.121015126908468, -0.292156360634725, 0.292156360634725, -0.121015126908468, -0.121015126908468, 0.292156360634725, -0.292156360634725, 0.121015126908468, 0.121015126908468, -0.292156360634725, 0.292156360634725, -0.121015126908467, -0.121015126908469, 0.292156360634725, -0.292156360634724, 0.121015126908467, 0.121015126908468, -0.292156360634725, 0.292156360634725, -0.121015126908468},
	{0.0977197537924274, -0.255833636800846, 0.316227766016838, -0.255833636800846, 0.0977197537924273, 0.0977197537924281, -0.255833636800847, 0.316227766016838, -0.255833636800846, 0.0977197537924271, 0.0977197537924277, -0.255833636800847, 0.316227766016838, -0.255833636800846, 0.0977197537924270, 0.0977197537924279, -0.255833636800847, 0.316227766016838, -0.255833636800846, 0.0977197537924268},
	{0.0738219058991409, -0.205373505467450, 0.292156360634725, -0.315252941349890, 0.269628494389924, -0.165228553867128, 0.0248109445657181, 0.121015126908469, -0.240461479767574, 0.307490367669328, -0.307490367669328, 0.240461479767574, -0.121015126908468, -0.0248109445657177, 0.165228553867129, -0.269628494389925, 0.315252941349890, -0.292156360634724, 0.205373505467449, -0.0738219058991401},
	{0.0494689214077114, -0.143564401525505, 0.223606797749979, -0.281761002650515, 0.312334477467278, -0.312334477467278, 0.281761002650515, -0.223606797749979, 0.143564401525505, -0.0494689214077107, -0.0494689214077114, 0.143564401525505, -0.223606797749979, 0.281761002650515, -0.312334477467278, 0.312334477467278, -0.281761002650515, 0.223606797749978, -0.143564401525503, 0.0494689214077104},
	{0.0248109445657177, -0.0738219058991411, 0.121015126908468, -0.165228553867129, 0.205373505467450, -0.240461479767575, 0.269628494389924, -0.292156360634725, 0.307490367669328, -0.315252941349890, 0.315252941349890, -0.307490367669328, 0.292156360634725, -0.269628494389923, 0.240461479767575, -0.205373505467449, 0.165228553867128, -0.121015126908468, 0.0738219058991390, -0.0248109445657143}
};

const std::vector<std::vector<float>> basisVectorsFloat =
{
	{0.223606797749979f, 0.223606797749979f, 0.223606797749979f, 0.223606797749979f, 0.223606797749979f, 0.223606797749979f, 0.223606797749979f, 0.223606797749979f, 0.223606797749979f, 0.223606797749979f, 0.223606797749979f, 0.223606797749979f, 0.223606797749979f, 0.223606797749979f, 0.223606797749979f, 0.223606797749979f, 0.223606797749979f, 0.223606797749979f, 0.223606797749979f, 0.223606797749979f},
	{0.315252941349890f, 0.307490367669328f, 0.292156360634725f, 0.269628494389924f, 0.240461479767574f, 0.205373505467449f, 0.165228553867129f, 0.121015126908468f, 0.0738219058991409f, 0.0248109445657177f, -0.0248109445657176f, -0.0738219058991408f, -0.121015126908468f, -0.165228553867129f, -0.205373505467449f, -0.240461479767574f, -0.269628494389924f, -0.292156360634725f, -0.307490367669328f, -0.315252941349890f},
	{0.312334477467278f, 0.281761002650515f, 0.223606797749979f, 0.143564401525505f, 0.0494689214077114f, -0.0494689214077113f, -0.143564401525505f, -0.223606797749979f, -0.281761002650515f, -0.312334477467278f, -0.312334477467278f, -0.281761002650515f, -0.223606797749979f, -0.143564401525505f, -0.0494689214077114f, 0.0494689214077113f, 0.143564401525505f, 0.223606797749979f, 0.281761002650515f, 0.312334477467278f},
	{0.307490367669328f, 0.240461479767574f, 0.121015126908468f, -0.0248109445657176f, -0.165228553867129f, -0.269628494389924f, -0.315252941349890f, -0.292156360634725f, -0.205373505467449f, -0.0738219058991411f, 0.0738219058991409f, 0.205373505467449f, 0.292156360634725f, 0.315252941349890f, 0.269628494389924f, 0.165228553867129f, 0.0248109445657176f, -0.121015126908468f, -0.240461479767574f, -0.307490367669328f},
	{0.300750477503773f, 0.185874017230092f, 1.93633660727019e-17f, -0.185874017230092f, -0.300750477503773f, -0.300750477503773f, -0.185874017230092f, -5.80900982181058e-17f, 0.185874017230092f, 0.300750477503773f, 0.300750477503773f, 0.185874017230092f, 9.68168303635097e-17f, -0.185874017230092f, -0.300750477503773f, -0.300750477503773f, -0.185874017230092f, -1.35543562508914e-16f, 0.185874017230092f, 0.300750477503773f},
	{0.292156360634725f, 0.121015126908468f, -0.121015126908468f, -0.292156360634725f, -0.292156360634725f, -0.121015126908468f, 0.121015126908468f, 0.292156360634725f, 0.292156360634725f, 0.121015126908468f, -0.121015126908468f, -0.292156360634725f, -0.292156360634725f, -0.121015126908468f, 0.121015126908468f, 0.292156360634725f, 0.292156360634725f, 0.121015126908468f, -0.121015126908467f, -0.292156360634725f},
	{0.281761002650515f, 0.0494689214077114f, -0.223606797749979f, -0.312334477467278f, -0.143564401525505f, 0.143564401525504f, 0.312334477467278f, 0.223606797749979f, -0.0494689214077115f, -0.281761002650515f, -0.281761002650515f, -0.0494689214077112f, 0.223606797749979f, 0.312334477467278f, 0.143564401525505f, -0.143564401525504f, -0.312334477467278f, -0.223606797749979f, 0.0494689214077109f, 0.281761002650515f},
	{0.269628494389924f, -0.0248109445657176f, -0.292156360634725f, -0.240461479767574f, 0.0738219058991409f, 0.307490367669328f, 0.205373505467449f, -0.121015126908468f, -0.315252941349890f, -0.165228553867129f, 0.165228553867129f, 0.315252941349890f, 0.121015126908468f, -0.205373505467449f, -0.307490367669328f, -0.0738219058991409f, 0.240461479767574f, 0.292156360634725f, 0.0248109445657175f, -0.269628494389924f},
	{0.255833636800846f, -0.0977197537924274f, -0.316227766016838f, -0.0977197537924275f, 0.255833636800846f, 0.255833636800847f, -0.0977197537924273f, -0.316227766016838f, -0.0977197537924275f, 0.255833636800846f, 0.255833636800847f, -0.0977197537924272f, -0.316227766016838f, -0.0977197537924276f, 0.255833636800846f, 0.255833636800847f, -0.0977197537924272f, -0.316227766016838f, -0.0977197537924277f, 0.255833636800846f},
	{0.240461479767574f, -0.165228553867129f, -0.292156360634725f, 0.0738219058991409f, 0.315252941349890f, 0.0248109445657179f, -0.307490367669328f, -0.121015126908468f, 0.269628494389924f, 0.205373505467450f, -0.205373505467449f, -0.269628494389924f, 0.121015126908468f, 0.307490367669328f, -0.0248109445657181f, -0.315252941349890f, -0.0738219058991410f, 0.292156360634724f, 0.165228553867129f, -0.240461479767574f},
	{0.223606797749979f, -0.223606797749979f, -0.223606797749979f, 0.223606797749979f, 0.223606797749979f, -0.223606797749979f, -0.223606797749979f, 0.223606797749979f, 0.223606797749979f, -0.223606797749979f, -0.223606797749979f, 0.223606797749979f, 0.223606797749979f, -0.223606797749979f, -0.223606797749979f, 0.223606797749979f, 0.223606797749978f, -0.223606797749979f, -0.223606797749980f, 0.223606797749979f},
	{0.205373505467449f, -0.269628494389924f, -0.121015126908468f, 0.307490367669328f, 0.0248109445657176f, -0.315252941349890f, 0.0738219058991411f, 0.292156360634725f, -0.165228553867129f, -0.240461479767575f, 0.240461479767574f, 0.165228553867129f, -0.292156360634725f, -0.0738219058991421f, 0.315252941349890f, -0.0248109445657169f, -0.307490367669328f, 0.121015126908467f, 0.269628494389924f, -0.205373505467448f},
	{0.185874017230092f, -0.300750477503773f, -5.80900982181058e-17f, 0.300750477503773f, -0.185874017230092f, -0.185874017230093f, 0.300750477503773f, 7.36003649626590e-16f, -0.300750477503773f, 0.185874017230092f, 0.185874017230092f, -0.300750477503773f, -8.52183846062801e-16f, 0.300750477503773f, -0.185874017230092f, -0.185874017230093f, 0.300750477503773f, -1.55102667445532e-16f, -0.300750477503773f, 0.185874017230092f},
	{0.165228553867129f, -0.315252941349890f, 0.121015126908468f, 0.205373505467450f, -0.307490367669328f, 0.0738219058991406f, 0.240461479767574f, -0.292156360634725f, 0.0248109445657181f, 0.269628494389924f, -0.269628494389924f, -0.0248109445657175f, 0.292156360634725f, -0.240461479767574f, -0.0738219058991411f, 0.307490367669328f, -0.205373505467450f, -0.121015126908469f, 0.315252941349890f, -0.165228553867129f},
	{0.143564401525505f, -0.312334477467278f, 0.223606797749979f, 0.0494689214077117f, -0.281761002650515f, 0.281761002650515f, -0.0494689214077120f, -0.223606797749980f, 0.312334477467278f, -0.143564401525504f, -0.143564401525504f, 0.312334477467278f, -0.223606797749979f, -0.0494689214077125f, 0.281761002650516f, -0.281761002650515f, 0.0494689214077106f, 0.223606797749980f, -0.312334477467278f, 0.143564401525503f},
	{0.121015126908468f, -0.292156360634725f, 0.292156360634725f, -0.121015126908468f, -0.121015126908468f, 0.292156360634725f, -0.292156360634725f, 0.121015126908468f, 0.121015126908468f, -0.292156360634725f, 0.292156360634725f, -0.121015126908467f, -0.121015126908469f, 0.292156360634725f, -0.292156360634724f, 0.121015126908467f, 0.121015126908468f, -0.292156360634725f, 0.292156360634725f, -0.121015126908468f},
	{0.0977197537924274f, -0.255833636800846f, 0.316227766016838f, -0.255833636800846f, 0.0977197537924273f, 0.0977197537924281f, -0.255833636800847f, 0.316227766016838f, -0.255833636800846f, 0.0977197537924271f, 0.0977197537924277f, -0.255833636800847f, 0.316227766016838f, -0.255833636800846f, 0.0977197537924270f, 0.0977197537924279f, -0.255833636800847f, 0.316227766016838f, -0.255833636800846f, 0.0977197537924268f},
	{0.0738219058991409f, -0.205373505467450f, 0.292156360634725f, -0.315252941349890f, 0.269628494389924f, -0.165228553867128f, 0.0248109445657181f, 0.121015126908469f, -0.240461479767574f, 0.307490367669328f, -0.307490367669328f, 0.240461479767574f, -0.121015126908468f, -0.0248109445657177f, 0.165228553867129f, -0.269628494389925f, 0.315252941349890f, -0.292156360634724f, 0.205373505467449f, -0.0738219058991401f},
	{0.0494689214077114f, -0.143564401525505f, 0.223606797749979f, -0.281761002650515f, 0.312334477467278f, -0.312334477467278f, 0.281761002650515f, -0.223606797749979f, 0.143564401525505f, -0.0494689214077107f, -0.0494689214077114f, 0.143564401525505f, -0.223606797749979f, 0.281761002650515f, -0.312334477467278f, 0.312334477467278f, -0.281761002650515f, 0.223606797749978f, -0.143564401525503f, 0.0494689214077104f},
	{0.0248109445657177f, -0.0738219058991411f, 0.121015126908468f, -0.165228553867129f, 0.205373505467450f, -0.240461479767575f, 0.269628494389924f, -0.292156360634725f, 0.307490367669328f, -0.315252941349890f, 0.315252941349890f, -0.307490367669328f, 0.292156360634725f, -0.269628494389923f, 0.240461479767575f, -0.205373505467449f, 0.165228553867128f, -0.121015126908468f, 0.0738219058991390f, -0.0248109445657143f}
};