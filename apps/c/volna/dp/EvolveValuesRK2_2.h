inline void EvolveValuesRK2_2(const double *dT, double *outConservative, //OP_RW, discard
            double *inConservative, //OP_READ, discard
            double *midPointConservative, //OP_READ, discard
            double *out) //OP_WRITE

{
  outConservative[0] = 0.5*(outConservative[0] * *dT + midPointConservative[0] + inConservative[0]);
  outConservative[1] = 0.5*(outConservative[1] * *dT + midPointConservative[1] + inConservative[1]);
  outConservative[2] = 0.5*(outConservative[2] * *dT + midPointConservative[2] + inConservative[2]);

  outConservative[0] = outConservative[0] <= EPS ? EPS : outConservative[0];
  outConservative[3] = inConservative[3];

  //call to ToPhysicalVariables inlined
  double TruncatedH = outConservative[0] < EPS ? EPS : outConservative[0];
  out[0] = outConservative[0];
  out[1] = outConservative[1] / TruncatedH;
  out[2] = outConservative[2] / TruncatedH;
  out[3] = outConservative[3];
}
