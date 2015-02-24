#!/home/osaka/.rvm/rubies/ruby-1.9.3-p551/bin/ruby

require 'json'

def norm_rand(mu=0,sigma=1.0)
  r1,r2=Random.rand(),Random.rand();
  Math::sqrt(-1*Math::log(r1))*Math::cos(2*Math::PI*r2)*sigma+mu
end

# x0=ARGV[0].to_f
# x1=ARGV[1].to_f
# x2=ARGV[2].to_f
# x3=ARGV[3].to_f
# x4=ARGV[4].to_f
# x5=ARGV[5].to_f
# x6=ARGV[6].to_f

parsed = JSON.load(open('./_input.json'))
x0 = parsed["x0"]
x1 = parsed["x1"]
x2 = parsed["x2"]
x3 = parsed["x3"]
x4 = parsed["x4"]
x5 = parsed["x5"]
x6 = parsed["x6"]
x7 = parsed["x7"]
x8 = parsed["x8"]
x9 = parsed["x9"]
x10 = parsed["x10"]
x11 = parsed["x11"]
x12 = parsed["x12"]
x13 = parsed["x13"]
x14 = parsed["x14"]

result = { "y" => 100*x0*100*x1+100*x3*x4+norm_rand(mu=0,sigma=0.3*x5+1)+100*x10+100*x0*x12**2*x13+x8*x9*x10*x11+100*x2+0.01*x14**2
 }

JSON.dump(result, open('./_output.json','w'))
