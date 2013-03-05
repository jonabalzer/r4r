#include"main.h"

int main(int argc, char *argv[]) {

    if(argc==1) {

        cout << "Type \"lmtest 0\" for LM-method on Rosebrock's function or \"lmtest 1\" for robust rew-eighted LM-method on Osbourne's function." << endl;
        return 1;

    }

    char* option = argv[1];

    if((unsigned int)option[0]==48) {  // LM demo on Rosenbrock function

        mat M(0,0);
        CPreconditioner<mat,vec,double> precond = CPreconditioner<mat,vec,double>(M);
        CIterativeSolver<mat,vec,double> solver = CIterativeSolver<mat,vec,double>(precond,2,1e-12,true);

        CRosenbrockFunction problem = CRosenbrockFunction();
        CLevenbergMarquardt<mat> lms(problem,solver,0);

        // access solution vector
        vec& model = problem.Get();
        model(0) = -2.1;
        model(1) = 1;

        // iterate
        lms.Iterate(100,1e-10,1e-10,false);

        cout << "Result:" << endl;
        cout << model << endl;

        cout << "Ground truth:" << endl;
        cout << problem.GetGroundTruth() << endl;

        cout << "Error:" << endl;
        cout << model - problem.GetGroundTruth() << endl;

    }
    else if((unsigned int)option[0]==49) { // reweighted LS on Osbourne function

        mat M(0,0);
        CPreconditioner<mat,vec,double> precond = CPreconditioner<mat,vec,double>(M);
        CIterativeSolver<mat,vec,double> solver = CIterativeSolver<mat,vec,double>(precond,5,1e-12,true);

        COsbourneFunction problem = COsbourneFunction();
        problem.DisturbSamplePoints(10,1);

        CLevenbergMarquardt<mat> lms(problem,solver,1e-8);

        // init
        vec& model = problem.Get();
        model(0) = 0.5;
        model(1) = 1.5;
        model(2) = -1;
        model(3) = 0.01;
        model(4) = 0.02;

        lms.Iterate(500,4,1e-15,false,true);

        cout << "Error RWLS:" << endl;
        cout << model - problem.GetGroundTruth() << endl;

        // compare with standard LM
        vec& weights = problem.GetWeights();
        weights.Ones();

        model(0) = 0.5;
        model(1) = 1.5;
        model(2) = -1;
        model(3) = 0.01;
        model(4) = 0.02;

        lms.Iterate(2000,1e-10,1e-10,true);

        cout << "Error LS:" << endl;
        cout << model - problem.GetGroundTruth() << endl;

    }
    else
        cout << "Invalid option." << endl;

    return 0;

}




CRosenbrockFunction::CRosenbrockFunction():
    CLeastSquaresProblem(2,2) {

}

void CRosenbrockFunction::ComputeResidual(vec& r) {

    r(0) = m_weights(0)*10*(m_model(1) - m_model(0)*m_model(0));
    r(1) = m_weights(1)*(1 - m_model(0));

}

void CRosenbrockFunction::ComputeResidualAndJacobian(vec& r, mat& J) {

    r(0) = m_weights(0)*(10*(m_model(1) - m_model(0)*m_model(0)));
    r(1) = m_weights(1)*(1 - m_model(0));

    J(0,0) = -m_weights(0)*20*m_model(0);
    J(0,1) = m_weights(0)*10;
    J(1,0) = -m_weights(1)*1;
    J(1,1) = m_weights(1)*0;

}

double COsbourneFunction::m_y[33] = {};
double COsbourneFunction::m_t[33] = {};

COsbourneFunction::COsbourneFunction():
    CLeastSquaresProblem(33,5) {

    for(size_t i=0; i<33; i++) {

        m_t[i] = 10*i;
        m_y[i] = EvalFunction(m_t[i]);

    }

}

void COsbourneFunction::ComputeResidual(vec& r) {

    for(size_t i=0; i<GetNumberOfDataPoints(); i++)
        r(i) = m_weights(i)*(m_y[i] - (m_model(0) + m_model(1)*exp(-m_model(3)*m_t[i]) + m_model(2)*exp(-m_model(4)*m_t[i])));

}

void COsbourneFunction::ComputeResidualAndJacobian(vec& r, mat& J) {

    for(size_t i=0; i<GetNumberOfDataPoints(); i++) {

        r(i) = m_weights(i)*((m_model(0) + m_model(1)*exp(-m_model(3)*m_t[i]) + m_model(2)*exp(-m_model(4)*m_t[i])) - m_y[i]);

        J(i,0) = m_weights(i)*1;
        J(i,1) = m_weights(i)*exp(-m_model(3)*m_t[i]);
        J(i,2) = m_weights(i)*exp(-m_model(4)*m_t[i]);
        J(i,3) = -m_weights(i)*m_t[i]*m_model(1)*exp(-m_model(3)*m_t[i]);
        J(i,4) = -m_weights(i)*m_t[i]*m_model(2)*exp(-m_model(4)*m_t[i]);

    }

}

void COsbourneFunction::DisturbSamplePoints(size_t noutlier, double strength) {

    srand(time(NULL));

    for(size_t i=0; i<noutlier; i++) {

        size_t index = (size_t)(((double)rand()/(double)RAND_MAX)*32);
        m_y[index] = m_y[index] + ((double)rand()/(double)RAND_MAX)*strength;

    }

}

void COsbourneFunction::PrintData() {

    for(size_t i=0; i<GetNumberOfDataPoints(); i++)
        cout << m_t[i] << " " << m_y[i] << endl;

}

double COsbourneFunction::EvalFunction(double t) {

    return 0.37541 + 1.93585*exp(-0.01287*t) - 1.46469*exp(-0.02212*t);

}

vec COsbourneFunction::GetGroundTruth() {

    vec x(5);
    x(0) = 0.37541;
    x(1) = 1.93585;
    x(2) = -1.46469;
    x(3) = 0.01287;
    x(4) = 0.02212;

    return x;

}

bool COsbourneFunction::SaveData(const char* filename) {

    mat data(33,2);

    for(size_t i=0; i<GetNumberOfDataPoints(); i++) {

        data(i,0) = m_t[i];
        data(i,1) = m_y[i];

    }

    ofstream out(filename);

    if(!out.is_open())
        return 1;

    out << data;

    out.close();

    return 0;

}


