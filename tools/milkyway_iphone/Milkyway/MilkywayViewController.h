#import <UIKit/UIKit.h>
#import <MessageUI/MessageUI.h>

extern double _milkywaySeparationGlobalProgress;
int _iphone_main(int argc, const char ** argv);

@interface MilkywayViewController : UIViewController<MFMailComposeViewControllerDelegate>
{
    NSTimer * progressTimer;
    NSString * logPath;
    NSThread * milkywayThread;
}

@property (nonatomic, retain) IBOutlet UIProgressView * progressView;
@property (nonatomic, retain) IBOutlet UITextView * outputView;
@property (nonatomic, retain) IBOutlet UILabel * statusLabel;
@property (nonatomic, retain) IBOutlet UIButton * submitButton;
@property (nonatomic, retain) IBOutlet UIButton * startButton;

- (NSString *)platform;

- (IBAction)start:(id)sender;
- (IBAction)submit:(id)sender;

@end
