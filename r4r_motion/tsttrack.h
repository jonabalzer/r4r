/*
 * tsthash.h
 *
 *  Created on: Mar 17, 2012
 *      Author: jbalzer
 */

#ifndef STTRACK_H_
#define STTRACK_H_


#define COMPUTE_ID_TST 0

#include "tracker.h"
#include "features.h"
//#include <tools/stree.h>

namespace R4R {


/*! \brief ST tracker based on combinatorial feature matching
 *
 * \details WARNING: Under development. May not work properly.
 *
 */
class CTST:public CTracker {

public:


	//! Constructor.
	CTST(CParameters params);

	//! Initializes the tree with the given image representing the root node.
	virtual bool Init(std::vector<cv::Mat>& pyramid);

	//! \copydoc CTracker::Update(cv::Mat&,cv::Mat&)
	virtual bool Update(std::vector<cv::Mat>& pyramid0, std::vector<cv::Mat>& pyramid1);

	//! \copydoc CTracker::Clean()
    virtual void Clean(std::vector<cv::Mat>& pyramid0, std::vector<cv::Mat>& pyramid1) { PruneTree(m_root); };

	//! \copydoc CTracker::AddTracklets(cv::Mat&)
	virtual bool AddTracklets(std::vector<cv::Mat>& pyramid);

	//! \copydoc CTracker::Draw(cv::Mat&)
	virtual void Draw(cv::Mat& img) { DrawChildren(img,m_root); };

	//! Saves the current tree to a file in the dot format.
	bool SaveTreeToFile(const char* filename);

protected:

	shared_ptr<CSTTracklet> m_root;														//!< root node of the tree
	cv::GoodFeaturesToTrackDetector m_detector;											//!< feature detector
	//cv::FastFeatureDetector m_detector;

	//! Writes the children of a node to a file stream.
	void WriteChildren(std::ofstream& os, shared_ptr<CSTTracklet> node);

	//! Draws the tree starting from a given node.
	void DrawChildren(cv::Mat& img, shared_ptr<CSTTracklet> node);

	//! Finds parent tracklet of a point.
	shared_ptr<CSTTracklet> FindParent(shared_ptr<CTracklet> tracklet);

	//! Constructs tree from latest state of all active tracklets.
	void ConstructTree();

	//! Removes all nodes below a certain node from the tree.
	bool RemoveSubTree(shared_ptr<CSTTracklet> node);

	//! Clears the entire tree.
	void ClearTree() { RemoveSubTree(m_root); };

	//! Remove inactive tracklets from the tree.
	void PruneTree(shared_ptr<CSTTracklet> node);

	//! Integrates orphan tracklets into the tree.
	void UpdateTree();

	//! Adds features to the trackers that adhere to distance condition.
	void Detect(std::vector<cv::Mat>& pyramid);

	//! Advances all features below a given node.
	void TrackSubTree(std::vector<cv::Mat>& pyramid0, std::vector<cv::Mat>& pyramid1, shared_ptr<CSTTracklet> node, vec t0);


};



/*! \brief node in a selection tree
 *
 *
 *
 */
class CTSTNode {

	friend class CTSTLK;
	friend class CTSTBTDTracker;

public:

	//! Standard constructor.
	CTSTNode();

	//! Constructor.
	CTSTNode(shared_ptr<CTracklet> tracklet);

	//! Destructor.
	virtual ~CTSTNode();

	//! Returns true if the node is a leave.
	bool IsLeave() { return m_children.size() == 0; };

protected:

	shared_ptr<CTracklet> m_tracklet;			//!< pointer to a tracklet
	CTSTNode* m_parent;							//!< pointer to parent
	std::list<CTSTNode*> m_children;			//!< pointers to children


};

/*! \brief ST tracker based on gray value SSD matching
 *
 * \details Implementation of the feature tracking algorithm described in [Lee2011]. The following parameters have
 * to be supplied by means of an object of CParameters:
 *
 * - SCALE number of levels in the selection tree
 * - FEATURE_THRESHOLD detection threshold
 * - TERMINATION_TIME number of frames to process in video stream
 * - REFRESH_RATE number of frames between redetections
 * - TRACKING_HSIZE half window size (for low-level motion estimation)
 * - MINIMAL_FEATURE_HDISTANCE half of the minimal distance that new features keep to existing ones
 * - LK_PYRAMID_LEVEL number of scales to use in low-level LK motion estimation
 * - MAX_ITER maximum number of Gauss-Newton steps to execute in LK motion estimation
 * - ACCURACY time derivative of objective below which to break off Gauss-Newton iteration
 * - LAMBDA relative weight of the spatial image derivatives impact to the optical flow estimation
 *
 * \todo Replace OpenCV rectangle with R4R CRectangle class. Eliminate external tree class, use specialized tracklet class instead.
 *
 */
class CTSTLK:public CTracker {

public:


	//! Constructor.
	CTSTLK(CParameters params);

	//! Destructor.
	~CTSTLK() { delete m_root; };

	//! Initializes the tree with the given image representing the root node.
	virtual bool Init(std::vector<cv::Mat>& pyramid);

	/*! \copydoc CTracker::Update(cv::Mat&,cv::Mat&)
	 *
	 * \details All nodes in the tree are collected by CTSTLK::CollectSiblings() independent of
	 * whether their status is valid or not. To get rid of bad tracks, invoke CTSTLK::Clean(),
	 * which effectively prunes the tree.
	 *
	 */
	virtual bool Update(std::vector<cv::Mat>& pyramid0, std::vector<cv::Mat>& pyramid1);

	/*! \copydoc CTracker::Clean()
	 *
	 *  \details First, CTSTLK::ValidateDescendants() is called to mark the nodes that violate
	 *  the membership condition of the tree. CTSTLK::Prune() removes all the nodes plus
	 *  the bad ones detected in CTSTLK::Update() and all of their descendants from the
	 *  tree.
	 *
	 */
    virtual void Clean(std::vector<cv::Mat>& pyramid0, std::vector<cv::Mat>& pyramid1);

	//! \copydoc CTracker::AddTracklets(cv::Mat&)
	virtual bool AddTracklets(std::vector<cv::Mat>& pyramid);

	//! \copydoc CTracker::Draw(cv::Mat&)
	//virtual void Draw(cv::Mat& img);

	//! Saves the current tree to a file in the dot format.
	bool SaveToFile(const char* filename);

protected:

	CTSTNode* m_root;																	//!< root node of the tree
//	cv::GoodFeaturesToTrackDetector m_detector;											//!< feature detector
	cv::FastFeatureDetector m_detector;

	//! Groups descendants of a node into groups of siblings.
	void CollectSiblings(std::vector<std::list<CTSTNode*> >& nodes, CTSTNode* node);

	//! Expands an existing tree starting from a give node.
	bool Reproduce(std::vector<cv::Mat>& pyramid, CTSTNode* node);

	/*! \brief Checks the descendants of a node whether they still fulfill the membership condition for the tree.
	 *
	 * \details Not all of the tree is traversed. If a node is found that violates the condition, the recursion
	 * is broken because all descendants will be unreliable and deleted when calling CTSTLK::Clean().
	 *
	 * \param[in] node node from where recursive evaluation starts
	 * \param[in] img current frame
	 *
	 */
	bool ValidateDescendants(cv::Mat& img, CTSTNode* node);

	//! Removes the node of a bad track and all of its descendants.
	bool Prune(CTSTNode* node);

};





}

#endif /* STTRACK_H_ */
