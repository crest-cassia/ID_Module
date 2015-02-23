#!/home/osaka/.rvm/rubies/ruby-1.9.3-p551/bin/ruby

require 'json'

parsed = JSON.load(open('./_input.json'))

#x0=ARGV[0].to_f
#x1=ARGV[1].to_f
#x2=ARGV[2].to_f
#x3=ARGV[3].to_f

x0 = parsed["x0"]
x1 = parsed["x1"]
x2 = parsed["x2"]
x3 = parsed["x3"]

def norm_rand(mu=0,sigma=1.0)
  r1,r2=Random.rand(),Random.rand();
  Math::sqrt(-1*Math::log(r1))*Math::cos(2*Math::PI*r2)*sigma+mu
end

result = { "y" => 10*x0*x1+norm_rand(mu=0,sigma=0.3*x2+1) }

JSON.dump(result, open('./_output.json','w'))
