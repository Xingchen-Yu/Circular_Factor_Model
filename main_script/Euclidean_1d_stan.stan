data {
  int<lower=0> I;   // number of items
  int<lower=0> J;   // number of responses
  matrix[I, J] y;   // predictor matrix
}
parameters {
  vector[I] beta;
  vector[J] alpha;
  vector[J] mu;
}
model {
  alpha ~ std_normal();         
  beta ~ std_normal();
  mu ~ std_normal();
  for(i in 1:I)
  y[i] ~  bernoulli(Phi(alpha * beta[i] + mu));
}

