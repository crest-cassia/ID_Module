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
x15 = parsed["x15"]
x16 = parsed["x16"]
x17 = parsed["x17"]
x18 = parsed["x18"]
x19 = parsed["x19"]

x20 = parsed["x20"]
x21 = parsed["x21"]
x22 = parsed["x22"]
x23 = parsed["x23"]
x24 = parsed["x24"]
x25 = parsed["x25"]
x26 = parsed["x26"]
x27 = parsed["x27"]
x28 = parsed["x28"]
x29 = parsed["x29"]

x30 = parsed["x30"]

     # Example 4 ver3

result = { "y" => 100*x7*100*x26+100*x4*100*x15+100*x18*100*x1+100*x4*100*x20+100*x11*100*x15+100*x0*100*x1+100*x4*100*x27+100*x7*100*x1+100*x21*100*x4+100*x14*100*x2+100*x24*100*x13+100*x7*100*x22+100*x24*100*x9+100*x4*100*x16+100*x23*100*x29+norm_rand(mu=0,sigma=1.0) }

JSON.dump(result, open('./_output.json','w'))
