#ifndef WORKERINCLUDED
#define WORKERINCLUDED
#include <iostream>
#include <semaphore>
#include <thread>

using namespace std;

// namespace fatrop
// {
//     template <typename WorkFunction>
//     class Worker
//     {
//         binary_semaphore signal2worker;
//         binary_semaphore signal2main;
//         // const WorkFunction wf;
//         WorkFunction wf;
//         atomic<bool> running;
//         thread t1;
//         void run()
//         {
//             while (true)
//             {
//                 signal2worker.acquire();
//                 if (running)
//                 {
//                     wf();
//                     signal2main.release();
//                 }
//                 else
//                 {
//                     signal2main.release();
//                     break;
//                 }
//             }
//         }

//     public:
//         Worker(WorkFunction &&wf) : signal2worker(0), signal2main(1), wf(wf), running(true), t1(&Worker::run, this)
//         {
//         }
//         /** \brief Wait for function evaluation to be completed */
//         void wait()
//         {
//             signal2main.acquire();
//             signal2main.release();
//         }
//         /** \brief Evaluate function */
//         void eval()
//         {
//             signal2main.acquire();
//             signal2worker.release();
//         }
//         /** \brief Finalize the worker, join threads */
//         void finalize()
//         {
//             signal2main.acquire();
//             running = false;
//             signal2worker.release();
//             t1.join();
//             signal2main.release();
//         }
//         ~Worker()
//         {
//             if (running == true)
//             {
//                 this->finalize();
//             }
//         }
//     };
// }
#endif // WORKERINCLUDED